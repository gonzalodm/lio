module lr_data
!Nvirt = number of virtual molecular orbitals
!dim = dimension of matrix Linear Response (Nvirt x NCO)
!cbas,cbasx = Needed for libint (Temporary)
!nstates = number of excited states to calculate
!root = excited state chosen for optimization
   implicit none

   logical :: lresp = .false.
   logical :: fbas = .true.
   integer :: Nvirt, dim, NCOlr, Mlr
   integer :: nstates = 3
   integer :: root = 0
   ! PASOS EN LR EN DINAMICA
   integer :: StepLR = 0
   real*8, dimension(:,:), allocatable :: eigenvec, cbas, cbasx
   real*8, dimension(:), allocatable :: eigenval
! Use Frozen Core Approximation
! nfo = number of molecular orbitals occupied with lowest energy deleted
! nfv = number of molecular orbiatls virtual with higher energy deleted
   logical :: FCA = .false.
   integer :: nfo = 3
   integer :: nfv = 3
! Matrix needed for change basis
   real*8, dimension(:,:), allocatable :: Coef_trans, Cocc, Cocc_trans, &
                                          Cvir, Cvir_trans
end module lr_data

module lrtddft
! This is a main subroutine for Linear Response calculate.
! This module performs a Tamm-Dancoff Aproximation (TDA)

contains

   subroutine linear_response(MatCoef,VecEne)
   use lr_data, only: Nvirt,dim,nstates,eigenval,&
                      eigenvec,cbas,root,FCA,nfo,nfv,&
                      NCOlr, Mlr
   use garcha_mod, only: NCO, M, c, a
   implicit none

   real*8, intent(in) :: MatCoef(M,M)
   real*8, intent(in) :: VecEne(M)

!  COUNTERS
   integer :: i, j, k
!  IF FCA IS USED
   real*8, dimension(:), allocatable :: Ene_LR
   real*8, dimension(:,:), allocatable :: Coef_LR
   real*8, dimension(:,:,:,:), allocatable :: K2eAO
!  DAVIDSON OPTIONS
   integer :: maxIter, iter ! ITERATION
   integer :: vec_dim, iv, first_vec, newvec ! VECTOR INDEXES
   integer :: max_subs, Subdim ! DIMENSION SUBSPACE
   real*8, dimension(:,:), allocatable :: tvecMO,vmatAO,Fv,AX,H,eigvec
   real*8, dimension(:,:), allocatable :: RitzVec,ResMat
   real*8, dimension(:), allocatable :: eigval, val_old, Osc
   real*8, dimension(:,:,:), allocatable :: tmatMO,tmatAO,FM2,FXC
   integer :: calc_2elec = 0
   logical :: conv = .false.

   call g2g_timer_start('LINEAR RESPONSE')

! INITIALIZATION OF VARIABLES IF FCA IS USED
   if (FCA .eqv. .true.) then
   print*,"Using Frozen Core Approximation"
   print*,"nfo, nfv", nfo, nfv
      Nvirt = M - NCO - nfv
      NCOlr = NCO - nfo
      Mlr = M - nfo - nfv
      dim = Nvirt * NCOlr
      allocate(Coef_LR(M,Mlr),Ene_LR(Mlr))
      do i=1, NCOlr
        Coef_LR(:,i) = MatCoef(:,i+nfo)
        Ene_LR(i) = VecEne(i+nfo)
      enddo
      do i=1, Nvirt
         Coef_LR(:,NCOlr+i) = MatCoef(:,i+NCO)
         Ene_LR(NCOlr+i) = VecEne(i+NCO)
      enddo
   else
      nfo = 0
      Mlr = M
      NCOlr = NCO
      Nvirt = M - NCO
      dim = Nvirt * NCO
      allocate(Coef_LR(M,M),Ene_LR(M))
      Coef_LR = MatCoef
      Ene_LR = VecEne
   endif

!  INITIALIZATION OF MATRIX NEEDED FOR CHANGE BASIS
   call basis_init(Coef_LR,M,NCO,Nvirt)

   print*, ""
   print*,"======================================="
   print*,"         LINEAR RESPONSE - TDA"
   print*,"======================================="
   print*, ""

!  CHECK INPUT VARIABLES
   if (nstates > dim) then
      print*, "NUMBER OF EXCITED STATES THAT YOU WANT IS BIGGER &
               THAN DIMENSION MATRIX"
      nstates = dim
      print*, "WE CALCULATED", nstates
   endif

!  DAVIDSON INITIALIZATION
   allocate(val_old(nstates),Osc(nstates))

   val_old = 1.0D0
   vec_dim = 4 * nstates
   first_vec = 0
   newvec = 0

   if(vec_dim >= dim) then
      vec_dim = dim
      Subdim = dim
      maxIter = 1
      max_subs = dim
      allocate(AX(dim,dim),tvecMO(dim,dim))
   else
      max_subs = 200
      if (max_subs > dim) max_subs = dim
      maxIter = 50
      Subdim = vec_dim
      allocate(AX(dim,max_subs),tvecMO(dim,max_subs))
   endif

!  INITIAL GUESS
   call vec_init(tvecMO,dim,vec_dim)

!  SAVE INTEGRALS COULOMP TYPE
   allocate(K2eAO(M,M,M,M),vmatAO(M,M),Fv(M,M))
   allocate(RitzVec(dim,nstates),ResMat(dim,nstates))
   K2eAO = 0.0D0; vmatAO = 0.0D0; Fv = 0.0D0; RitzVec = 0.0D0

!  PRINT INFORMATION
   write(*,"(1X,A,22X,I2,2X,I2,2X,I2)") "NCO, NVIRT, M",NCO,Nvirt,M
   write(*,"(1X,A,11X,I3)") "DIMENSION OF FULL MATRIX",dim
   write(*,"(1X,A,23X,I3)") "MAX SUBSPACE",max_subs
   write(*,"(1X,A,20X,I2)") "MAX ITERATIONES",maxIter
   write(*,"(1X,A,3X,I4)") "NUMBER OF INITIAL TRIALS VECTORS",vec_dim

!  DAVIDSON START
   do iter=1,maxIter
      write(*,"(A)") " "
      write(*,"(1X,A,6X,I2)") "ITERATION:",iter
      write(*,"(1X,A,7X,I4)") "SUBSPACE:",Subdim
      write(*,"(1X,A,1X,I4)") "VECTORS INSIDE:",vec_dim
      write(*,"(1X,A,4X,I2)") "LAST VECTOR:",first_vec

!    CHANGE VECTORS TO MATRIX IN OM
     if(allocated(tmatMO)) deallocate(tmatMO)
        allocate(tmatMO(M,M,vec_dim))
     call vecMOtomatMO(tvecMO,tmatMO,M,NCO,Nvirt, &
                       Subdim,vec_dim,first_vec,dim)

!    CHANGE BASIS MO -> AO FOR EACH VECTOR
     if(allocated(tmatAO)) deallocate(tmatAO)
        allocate(tmatAO(M,M,vec_dim))
     call matMOtomatAO(tmatMO,tmatAO,Coef_LR,M,vec_dim,.true.)

!    CALCULATE 2E INTEGRALS
     allocate(FM2(M,M,vec_dim)); FM2 = 0.0D0
     call g2g_calculate2e(tmatAO,K2eAO,cbas,vec_dim,FM2,calc_2elec)
     calc_2elec = 1 ! TURN OFF CALCULATE INTEGRALS

!    CALCULATE DFT INTEGRALS
     allocate(FXC(M,M,vec_dim))
     do iv=1,vec_dim ! LOOPS VECTORS INSIDE OF CYCLE
        call formred(tmatMO,Coef_LR,vmatAO,M,dim,vec_dim,iv) ! FORM T * C**T
        call g2g_calculatedft(vmatAO,Coef_LR,Fv,NCO)
        FXC(:,:,iv) = Fv
        Fv = 0.0D0
     enddo ! END LOOPS VECTORS

!    NOW FM2 CONTAIN THE 2e- AND DFT PARTS
     FM2 = FM2 + FXC 
     deallocate(FXC)

!    WE OBTAIN AX (NDIM x NVEC)
     call MtoIANV(FM2,Coef_LR,AX,M,NCO,dim, &
                  Subdim,vec_dim,first_vec)
     deallocate(FM2)

!    ADD (Ea - Ei)*Xia TO AX
     call addInt(AX,Ene_LR,tvecMO,dim,M, &
                 Subdim,NCO,vec_dim,first_vec)

!    WE OBTAIN SUBSPACE MATRIX
     allocate(H(Subdim,Subdim))
     call subspaceMat(AX,tvecMO,H,dim,Subdim)

!    DIAGONALIZATION OF SUBSPACE MATRIX
     if(allocated(eigvec)) deallocate(eigvec)
     if(allocated(eigval)) deallocate(eigval)
       allocate(eigvec(Subdim,Subdim),eigval(Subdim))
     call diagonH(H,Subdim,eigval,eigvec)
     deallocate(H)

!    CHANGE SUBSPACE TO FULL SPACE - OBTAIN RITZ VECTOR
     call RitzObtain(tvecMO,eigvec,RitzVec,dim,Subdim,nstates)

     if (maxIter == 1) then
       call OscStr(RitzVec,eigval,Coef_LR,Osc,M,NCO,Nvirt,dim,nstates)
       call PrintResults(RitzVec,eigval,Osc,dim,nstates)
       exit
     endif

!    OBTAIN RESIDUAL VECTORS
     call residual(eigval,eigvec,RitzVec,AX,ResMat, &
                   dim,Subdim,nstates)

!    CHECK CONVERGENCE AND ADD NEW VECTORS
     conv = .false. ; newvec = 0
     call new_vectors(ResMat,eigval,Ene_LR,tvecMO,val_old, &
                      dim,Subdim,nstates,M,NCO,newvec,conv)

     if(conv .eqv. .true.) then
       write(*,"(1X,A,I2,1X,A)") "CONVERGED IN:",iter,"ITERATIONS"
       call OscStr(RitzVec,eigval,Coef_LR,Osc,M,NCO,Nvirt,dim,nstates)
       call PrintResults(RitzVec,eigval,Osc,dim,nstates)
       exit
     else
!      ACTUALIZATION OF VECTORS INDEXES
       first_vec = Subdim 
       Subdim = Subdim + newvec
       vec_dim = newvec
     endif
   enddo !END DAVIDSON

   deallocate(vmatAO,Fv,ResMat,val_old,Osc,eigval, &
              eigvec,AX,tvecMO,tmatAO)

   if (root > 0 ) then
     write(*,"(1X,A,I2)") "FORM RELAXED DENSITY FOR EXICTED STATE:",root
     call Zvector(Coef_LR,Ene_LR,RitzVec(:,root),K2eAO,NCO,M,dim)
   endif
   
   deallocate(K2eAO,RitzVec,Coef_LR,Ene_LR)
   call basis_deinit()

   end subroutine linear_response

   subroutine vec_init(Vec,N,vecnum)
      implicit none

      integer, intent(in) :: N, vecnum
      real*8, intent(out) :: Vec(N,vecnum)

      integer :: i
  
      Vec = 0.0D0
      do i=1,vecnum
         Vec(i,i) = 1.0D0
      enddo
   end subroutine vec_init

   subroutine vecMOtomatMO(vec,mat,M,NCO,Nvirt,Sdim,Nvec,V1,dim)
      implicit none
  
      integer, intent(in) :: M, NCO, Nvirt, Nvec, V1, dim, Sdim
      real*8, intent(in) :: vec(dim,Sdim)
      real*8, intent(out) :: mat(M,M,Nvec)
  
      integer :: i,row,col,NCOc

      NCOc = NCO - 1
      mat = 0.0D0
     
      do i=1,Nvec
        do row=0,NCOc
        do col=1,Nvirt
           mat(NCOc-row+1,NCO+col,i) = vec(row*Nvirt+col,V1+i)
        enddo
        enddo
      enddo
   end subroutine vecMOtomatMO

   subroutine matMOtomatAO(MatMO,MatAO,C,M,numvec,transp)
   use lr_data, only: Coef_trans

      implicit none

      integer, intent(in) :: M, numvec
      logical, intent(in) :: transp
      real*8, intent(in) :: MatMO(M,M,numvec)
      real*8, intent(in) :: C(M,M)
      real*8, intent(out) :: MatAO(M,M,numvec)

      integer :: i, j, k
      real*8, dimension(:,:), allocatable :: scratch

!     CHANGE BASIS MO -> AO
      allocate(scratch(M,M))
      do i=1,numvec
         scratch = matmul(MatMO(:,:,i),Coef_trans)
         MatAO(:,:,i) = matmul(C,scratch)
      enddo

!     WE TRANSPOSE MATRIX FOR USE IT IN C
      if ( transp ) then
        do i=1,numvec
          scratch = MatAO(:,:,i)
          do j=1,M
          do k=1,M
            MatAO(k,j,i) = scratch(j,k)
          enddo
          enddo
        enddo
      endif

      deallocate(scratch)
   end subroutine matMOtomatAO

   subroutine formred(TmatMO,C,VecMOAO,M,dim,numvec,iv)
   use lr_data, only: Coef_trans

     implicit none

     integer, intent(in) :: M, dim, iv, numvec
     real*8, intent(in) :: TmatMO(M,M,numvec)
     real*8, intent(in) :: C(M,M)
     real*8, intent(out) :: VecMOAO(M,M)

     integer :: i, j

     vecMOAO = matmul(TmatMO(:,:,iv),Coef_trans(:,:))
   end subroutine formred

   subroutine MtoIANV(F,C,A,M,NCO,Ndim,Sdim,Nvec,V1)
   use lr_data, only: Cocc_trans

      implicit none

      integer, intent(in) :: M, Ndim, Nvec, NCO, V1, Sdim
      real*8, intent(in) :: F(M,M,Nvec), C(M,M)
      real*8, intent(inout) :: A(Ndim,Sdim)

      integer :: i, j, k, iv, row, Nvirt, NCOc
      real*8, dimension(:,:), allocatable :: B
 
      real*8 :: temp

      Nvirt = M - NCO
      NCOc = NCO + 1
      allocate(B(NCO,M))

      temp = 0.0D0
      do iv=1,Nvec
        B = matmul(Cocc_trans,F(:,:,iv))
        do i=1,NCO
        do j=NCOc,M
          do k=1,M
            temp = temp + B(NCOc-i,k) * C(k,j)
          enddo
          row = (i-1) * Nvirt + (j-NCO)
          A(row,V1+iv) = temp
          temp = 0.0D0
        enddo
        enddo
      enddo

      deallocate(B)
   end subroutine MtoIANV

   subroutine addInt(A,Ene,Vec,Ndim,M,Sdim,NCO,Nvec,V1)
      implicit none
      
      integer, intent(in) :: Ndim, M, NCO, Nvec, V1, Sdim
      real*8, intent(in) :: Ene(M), Vec(Ndim,Sdim)
      real*8, intent(inout) :: A(Ndim,Sdim)

      integer :: iv, i, j, Nvirt, NCOc, row

      Nvirt = M - NCO
      NCOc = NCO + 1
      do iv=1,Nvec
        do i=1,NCO
        do j=NCOc,M
          row = (i-1) * Nvirt + (j-NCO)
          A(row,V1+iv) = A(row,V1+iv) + (Ene(j) - Ene(NCOc-i)) * Vec(row,V1+iv)
        enddo
        enddo
      enddo
   end subroutine addInt

   subroutine subspaceMat(A,V,H,Ndim,Sdim)
      implicit none
 
      integer, intent(in) :: Ndim, Sdim
      real*8, intent(in) :: A(Ndim,Sdim), V(Ndim,Sdim)
      real*8, intent(out) :: H(Sdim,Sdim)

      integer :: i, j, k
      real*8, dimension(:,:), allocatable :: VT

      allocate(VT(Sdim,Ndim))
      do i=1,Ndim
      do j=1,Sdim
         VT(j,i) = V(i,j)
      enddo
      enddo

      H = matmul(VT,A)
      deallocate(VT)
   end subroutine subspaceMat

   subroutine diagonH(H,N,Val,Vec)
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: H(N,N)
   real*8, intent(out) :: Val(N),Vec(N,N)

   real*8, dimension(:), allocatable :: WORK
   integer, dimension(:), allocatable :: IWORK
   real*8 :: WI(N)
   integer :: i, j, info, LWORK, LIWORK

   Vec = H
   LWORK = -1
   LIWORK = -1

   allocate(WORK(1),IWORK(1))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   LWORK=WORK(1)
   LIWORK=IWORK(1)
   deallocate(WORK,IWORK); allocate(WORK(LWORK),IWORK(LIWORK))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   deallocate(WORK,IWORK)
   end subroutine diagonH

   subroutine PrintResults(vec,val,O,N,nstat)
   use lr_data, only: Mlr, NCOlr, nfo
   implicit none

   integer, intent(in) :: N, nstat
   real*8, intent(in) :: vec(N,nstat),val(nstat),O(nstat)

   integer :: i,j,from,to
   real*8 :: value_X

   from = NCOlr
   to = NCOlr + 1

   do j=1, nstat
   write(*,100) j, val(j), 45.56335D0/val(j), O(j)
   do i=1, N
      value_X = vec(i,j) / dsqrt(2.0D0)
      if ( abs(value_X) > 0.1D0 ) then
         write(*,101) from+nfo, to+nfo, value_X
      endif
      to = to + 1
      if ( to == Mlr+1 ) then
          from = from - 1
          to = NCOlr + 1
      endif
   enddo
      print*, " "
      from = NCOlr
      to = NCOlr + 1
   enddo

   100 FORMAT(1X,"STATE ",I2,3X,"ENERGY=",F8.4," Hartree, ",F12.6," nm"," OSC=",F8.4)
   101 FORMAT(6X,I2," -> ",I2,2X,F14.7)
   end subroutine PrintResults

   subroutine RitzObtain(V,Vsmall,R,Ndim,Sdim,Nstat)
      implicit none

      integer, intent(in) :: Ndim, Sdim, Nstat
      real*8, intent(in) :: V(Ndim,Sdim), Vsmall(Sdim,Nstat)
      real*8, intent(out) :: R(Ndim,Nstat)
   
      R = matmul(V,Vsmall)
   end subroutine RitzObtain

   subroutine residual(Val,Vec,R,A,W,Ndim,Sdim,Nstat)
      implicit none

      integer, intent(in) :: Ndim, Sdim, Nstat
      real*8, intent(in) :: Val(Nstat), Vec(Sdim,Nstat), R(Ndim,Nstat)
      real*8, intent(in) :: A(Ndim,Sdim)
      real*8, intent(out) :: W(Ndim,Nstat)

      integer :: i
      real*8, dimension(:), allocatable :: temp

      allocate(temp(Ndim))
      do i=1,Nstat
        temp = matmul(A(1:Ndim,1:Sdim),Vec(1:Sdim,i)) - Val(i) * R(1:Ndim,i)
        W(1:Ndim,i) = temp
      enddo
      deallocate(temp)
   end subroutine residual

   subroutine norma(V,N,norm)
      implicit none

      integer, intent(in) :: N
      real*8, intent(in) :: V(N)
      real*8,intent(out) :: norm

      integer :: i

      norm = 0.0d0
      do i=1,N
        norm = norm + V(i)*V(i)
      enddo
    end subroutine norma

    subroutine new_vectors(R,Esub,Eall,T,Eold,Ndim,Sdim,Nstat,M,NCO,New,conv)
       implicit none

       integer, intent(in) :: Ndim, Sdim, Nstat, M, NCO
       real*8, intent(in) :: R(Ndim,Nstat), Esub(Nstat), Eall(M)
       real*8, intent(inout) :: T(Ndim,Sdim+Nstat), Eold(Nstat)
       logical, intent(out) :: conv
       integer, intent(out) :: New

       integer :: i, iv, occ, virt, Nvirt, ind, NCOc
       real*8 :: temp, ERROR, MAX_ERROR, MAX_ENE, tolv, tole, diffE
       real*8, dimension(:,:), allocatable :: Qvec
       integer, dimension(:), allocatable :: valid_id

       allocate(Qvec(Ndim,Nstat))
       allocate(valid_id(Nstat))

       MAX_ERROR = 0.0D0
       MAX_ENE = 0.0D0
       tolv = 1.0D-6
       tole = 1.0D-6
       Nvirt = M - NCO
       NCOc = NCO + 1

       do iv=1,Nstat
         diffE = abs(Eold(iv) - Esub(iv))
         if(diffE > MAX_ENE) MAX_ENE = diffE
         write(*,"(1X,A,I2,1X,A,F14.7)") "Vector:",iv,"Change in energy:",diffE

         call norma(R(:,iv),Ndim,ERROR)
         if(ERROR > MAX_ERROR) MAX_ERROR = ERROR

         print*, ERROR
         if(ERROR > tolv .or. diffE > tole) then
            New = New + 1
            valid_id(New) = iv
            do i=1,Ndim
              ind = i - 1
              occ = NCO - (ind/Nvirt)
              virt = mod(ind,Nvirt) + NCOc
              temp = Eall(virt) - Eall(occ)
              temp = 1.0D0 / (Esub(iv) - temp)
              Qvec(i,iv) = temp * R(i,iv)
            enddo
         else
            write(*,"(1X,A,I2,1X,A)") "Vector:",iv,"Converged"
         endif
       enddo

       if(Sdim + New >= Ndim) then
          conv = .true.
          return
       endif

       print*,"VEC,ENE",MAX_ERROR,MAX_ENE
       if(MAX_ERROR < tolv .and. MAX_ENE < tole) then
         conv = .true.
       else
         conv = .false.
         do iv=1,New
           T(:,Sdim+iv) = Qvec(:,valid_id(iv))
         enddo
         call QRfactorization(T,Ndim,Sdim+New)
       endif

       Eold = Esub
       deallocate(Qvec,valid_id)
   end subroutine new_vectors

   subroutine QRfactorization(A,N,M)
      implicit none

      integer, intent(in) :: N, M
      real*8, intent(inout) :: A(N,M)

      real*8,dimension(:),allocatable :: WORK, TAU
      integer :: i,j,LWORK,INFO,K
      real*8 :: norma, ortog

      allocate(WORK(1),TAU(M))
      call dgeqrf(N,M,A,N,TAU,WORK,-1,INFO)
      LWORK = WORK(1)
      deallocate(WORK); allocate(WORK(LWORK))
      call dgeqrf(N,M,A,N,TAU,WORK,LWORK,INFO)
      call dorgqr(N,M,M,A,N,TAU,WORK,LWORK,INFO)
   end subroutine QRfactorization

   subroutine OscStr(X,Ene,Coef,OsSt,M,NCO,Nvirt,Ndim,Nstat)
      implicit none

      integer, intent(in) :: M, Ndim, Nstat, NCO, Nvirt
      real*8, intent(in) :: X(Ndim,Nstat),Ene(Nstat)
      real*8, intent(in) :: Coef(M,M)
      real*8, intent(out) :: OsSt(Nstat)

      integer :: i, j
      real*8, dimension(:,:), allocatable :: Tdip
      real*8, dimension(:,:,:), allocatable :: TdensAO,TdensMO

      allocate(Tdip(Nstat,3),TdensMO(M,M,Nstat),TdensAO(M,M,Nstat))

!     FORM TRANSITION DENSITY IN MO BASIS
      call vecMOtomatMO(X,TdensMO,M,NCO,Nvirt,Nstat,Nstat,0,Ndim)
!     CHANGE BASIS MO -> AO
      call matMOtomatAO(TdensMO,TdensAO,Coef,M,Nstat,.true.)
      deallocate(TdensMO)
!     CALCULATE TRANSITION DIPOLE
      call TransDipole(TdensAO,Tdip,M,Nstat)
      deallocate(TdensAO)
!     CALCULATE THE OSCILATOR STRENGHT
      call ObtainOsc(Tdip,Ene,OsSt,Nstat)
      deallocate(Tdip)
   end subroutine OscStr

   subroutine TransDipole(Tdens,Tdip,M,Nstat)
      use garcha_mod, only: RMM
      implicit none

      integer, intent(in) :: M, Nstat
      real*8, intent(inout) :: Tdens(M,M,Nstat)
      real*8, intent(out) :: Tdip(Nstat,3)

      integer :: i, j, ist
      real*8, dimension(:,:), allocatable :: RhoRmm

!     SAVE RHO FROM RMM
      allocate(RhoRmm(M,M))
      call spunpack_rho('L',M,RMM,RhoRmm)

      do ist=1,Nstat
         do i=1,M
         do j=1,i-1
            Tdens(i,j,ist) = Tdens(i,j,ist) + Tdens(j,i,ist)
         enddo
         enddo
         call sprepack('L',M,RMM,Tdens(:,:,ist))
         call dip(Tdip(ist,:))
      enddo
      Tdip = Tdip * 2.0D0 / dsqrt(2.0D0)

!     COPY RHO OLD INTO RMM
      do i=1,M
      do j=1,i-1
         RhoRMM(i,j) = RhoRMM(i,j) * 2.0D0
      enddo
      enddo
      call sprepack('L',M,RMM,RhoRMM)
      deallocate(RhoRMM)
   end subroutine TransDipole   
   
   subroutine ObtainOsc(dip,E,O,N)
      implicit none

      integer, intent(in) :: N
      real*8, intent(in) :: dip(N,3), E(N)
      real*8, intent(out) :: O(N)

      integer :: i, j
      real*8 :: dostres = 2.0D0 / 3.0D0
      real*8 :: temp = 0.0D0

      do i=1,N
      do j=1,3
         temp = temp + dip(i,j) * dip(i,j) * E(i) * dostres
      enddo
      O(i) = temp; temp = 0.0D0
      enddo
   end subroutine ObtainOsc

   subroutine Zvector(C,Ene,X,K4cen,NCO,M,Ndim)
   use lr_data, only: cbas
      implicit none

      integer, intent(in) :: NCO, M, Ndim
      real*8, intent(in) :: C(M,M), Ene(M)
      real*8, intent(in) :: X(Ndim), K4cen(M,M,M,M)

      integer :: i , j, Nvirt
      real*8, dimension(:,:), allocatable :: TundAO, Xmat, Gxc, TundMO
      real*8, dimension(:,:), allocatable :: FX, FT
      real*8, dimension(:,:,:), allocatable :: F2e, PA
      real*8, dimension(:,:), allocatable :: FXAB, FXIJ, FTIA, GXCIA
      real*8, dimension(:), allocatable :: Rvec

      Nvirt = M - NCO

!     FORM UNRELAXED DIFFERENCE DENSITY MATRIX
      allocate(TundMO(M,M))
      call UnDiffDens(X,TundMO,NCO,Nvirt,M,Ndim)

!     CHANGE BASIS T MO -> AO
      allocate(TundAO(M,M))
      call matMOtomatAO(TundMO,TundAO,C,M,1,.false.)
      deallocate(TundMO)

!     FORM TRANSITION DENSITY MATRIX
      allocate(Xmat(M,M)); Xmat = 0.0D0
      call XmatForm(X,C,Xmat,Ndim,NCO,Nvirt,M)

!     CALCULATE THIRD DERIVATIVE FOCK
      allocate(Gxc(M,M)); Gxc = 0.0D0
      call g2g_calculateg(Xmat,Gxc,3) 

!     CALCULATE SECOND DERIVATIVE FOCK
      allocate(FX(M,M)); FX = 0.0D0
      call g2g_calculateg(Xmat,FX,2) 
      allocate(FT(M,M)); FT = 0.0D0
      call g2g_calculateg(TundAO,FT,2) 

!     CALCULATE TWO ELECTRON FOCK
!     1 = tiene la parte 2e de Xmat
!     2 = tiene la parte 2e de Tund
      allocate(PA(M,M,2),F2e(M,M,2)); F2e = 0.0D0
      PA(:,:,1) = Xmat; PA(:,:,2) = TundAO
      deallocate(Xmat)
      call g2g_calculate2e(PA,K4cen,cbas,2,F2e,1)
      F2e = 2.0D0 * F2e

      ! FT + F2eT
      FT = F2e(:,:,2) + 2.0D0 * FT

      ! FX + F2eX
      FX = F2e(:,:,1) + 2.0D0 * FX
      deallocate(F2e)

!     CHANGE BASIS OF ALL FOCK TYPE MATRIX
      allocate(FXAB(Nvirt,Nvirt),FXIJ(NCO,NCO))
      allocate(FTIA(NCO,Nvirt),GXCIA(NCO,Nvirt))
      call ChangeBasisF(FX,FT,Gxc,C,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
      deallocate(FX,FT,Gxc)

!     CALCULATE VECTOR R (A * X = R)
      allocate(Rvec(Ndim))
      call RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,NCO,Nvirt,Ndim)
 
!     SOLVE EQUATION AX=R WITH PCG METHOD     
      call PCG_solve(Rvec,K4cen,TundAO,C,Ene,M,NCO,Nvirt,Ndim)
      
      deallocate(TundAO,FXAB,FXIJ,FTIA,GXCIA,Rvec)
   end subroutine Zvector

   subroutine UnDiffDens(X,T,NCO,Nvirt,M,Ndim)
      implicit none

      integer, intent(in) :: NCO, Nvirt, M, Ndim
      real*8, intent(in) :: X(Ndim)
      real*8, intent(out) :: T(M,M)

      integer :: i, j, pos, NCOc
      real*8 :: raiz2
      real*8, dimension(:,:), allocatable :: XM, XMtrans, Ptrash

      allocate(XM(NCO,Nvirt),XMtrans(Nvirt,NCO))

!     WE NORMALIZE TO 1/2
      raiz2 = 1.0D0 / dsqrt(2.0D0)
      T = 0.0D0
      NCOc = NCO + 1

      do i=1,NCO
      do j=1,Nvirt
        pos = (i-1) * Nvirt + j
        XM(i,j) = X(pos) * raiz2
        XMtrans(j,i) = X(pos) * raiz2
      enddo
      enddo

!     FORM UNRELAXED DIFFERENCE DENSITY MATRIX
      allocate(Ptrash(NCO,NCO))
      ! FORM BLOCK OCC-OCC
      Ptrash = (-1.0D0) * matmul(XM,XMtrans)
      do i=1,NCO
      do j=i,NCO
         T(NCOc-i,NCOc-j) = Ptrash(i,j)
         T(NCOc-j,NCOc-i) = Ptrash(i,j)
      enddo
      enddo

      ! FORM BLOCK VIR - VIR
      deallocate(Ptrash); allocate(Ptrash(Nvirt,Nvirt))
      Ptrash = matmul(XMtrans,XM)
      do i=1,Nvirt
      do j=i,Nvirt
         T(NCO+i,NCO+j) = Ptrash(i,j)
         T(NCO+j,NCO+i) = Ptrash(i,j)
      enddo
      enddo

      deallocate(Ptrash,XM,XMtrans)
   end subroutine UnDiffDens

   subroutine XmatForm(Vec,Coef,Mat,Ndim,NCO,Nvirt,M)
   use lr_data, only: Coef_trans
      implicit none

      integer, intent(in) :: Ndim, NCO, Nvirt, M
      real*8, intent(in) :: Vec(Ndim), Coef(M,M)
      real*8, intent(out) :: Mat(M,M)

      integer :: NCOc, row, col, pos
      real*8 :: raiz2
      real*8, dimension(:,:), allocatable ::SCR

      NCOc = NCO + 1
      raiz2 = 1.0D0 / dsqrt(2.0D0)

      do row=1,NCO
      do col=1,Nvirt
        pos = (row-1) * Nvirt + col
        Mat(NCOc-row,NCO+col) = Vec(pos) * raiz2
      enddo
      enddo

      allocate(SCR(M,M))
      SCR = matmul(Coef,Mat)
      Mat = matmul(SCR,Coef_trans)

      deallocate(SCR)
   end subroutine XmatForm

   subroutine ChangeBasisF(FX,FT,Gxc,C,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
   use lr_data, only: Cocc, Cocc_trans, Cvir, Cvir_trans

      implicit none

      integer, intent(in) :: M, Nvirt, NCO
      real*8, intent(in) :: C(M,M)
      real*8, intent(inout) :: FX(M,M), FT(M,M), Gxc(M,M)
      real*8, intent(out) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
      real*8, intent(out) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)

      integer :: i, j
      real*8, dimension(:,:), allocatable :: scratch

!     FORM FX IN BASIS VIRT X VIRT
      allocate(scratch(M,Nvirt))
      scratch = matmul(FX,Cvir)
      FXAB = matmul(Cvir_trans,scratch)
      deallocate(scratch)
 
!     FORM FX IN BASIS OCC X OCC
      allocate(scratch(M,NCO))
      scratch = matmul(FX,Cocc)
      FXIJ = matmul(Cocc_trans,scratch)
      deallocate(scratch)

!     FORM FT IN BASIS OCC X VIR
      allocate(scratch(M,Nvirt))
      scratch = matmul(FT,Cvir)
      FTIA = matmul(Cocc_trans,scratch)
      deallocate(scratch)

!     FORM GXC IN BASIS OCC X VIR
      allocate(scratch(M,Nvirt))
      scratch = matmul(Gxc,Cvir)
      GXCIA = matmul(Cocc_trans,scratch)
      deallocate(scratch)
   end subroutine ChangeBasisF

   subroutine RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,NCO,Nvirt,Ndim)
      implicit none
     
      integer, intent(in) :: NCO, Nvirt, Ndim
      real*8, intent(in) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
      real*8, intent(in) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)
      real*8, intent(in) :: X(Ndim)
      real*8, intent(out) :: Rvec(Ndim)

      integer :: i, a, b, j, posf, pos1, NCOc
      real*8 :: temp1, temp2, raiz2
 
      temp1 = 0.0D0; temp2 = 0.0D0
      raiz2 = 1.0D0 / dsqrt(2.0D0)
      NCOc = NCO + 1

      do i=1,NCO
      do a=1,Nvirt
         posf = (i-1) * Nvirt + a 
         ! VIRTUAL PART
         do b=1,Nvirt
            pos1 = (i-1) * Nvirt + b
            temp1 = temp1 + X(pos1) * FXAB(a,b) * raiz2
         enddo
         ! OCC PART
         do j=1,NCO
            pos1 = (j-1) * Nvirt + a 
            temp2 = temp2 + X(pos1) * FXIJ(NCOc-i,NCOc-j) * raiz2
         enddo
         ! OBRAIN RIA IN VECTOR FORM
         Rvec(posf) = temp2 - (temp1 + FTIA(NCOc-i,a) + 2.0D0*GXCIA(NCOc-i,a))
         temp1 = 0.0D0; temp2 = 0.0D0
      enddo
      enddo
   end subroutine RCalculate
 
   subroutine PCG_solve(bvec,K4,Rho_urel,Coef,E,M,NCO,Nvirt,Ndim)
      use lr_data, only: cbas
      implicit none

      integer, intent(in) :: NCO, Nvirt, Ndim, M
      real*8, intent(in) :: bvec(Ndim), E(M), Coef(M,M), Rho_urel(M,M)
      real*8, intent(in) :: K4(M,M,M,M)

      integer :: i, j, iter, maxIter
      real*8 :: beta, alpha
      logical :: conv = .false.
      real*8, dimension(:), allocatable :: R, Z, Pk, Mprec, ApIA, X
      real*8, dimension(:,:), allocatable :: Pmat, F2e, Fxc, Ftot

!     START PRECONDITINED CONJUGATE GRADIENT
      maxIter = 50

!     INITIAL GUESS: Xo = 0
      allocate(R(Ndim))
      R = bvec

!     CALCULATE PRECONDITIONED M^(-1)
      allocate(Mprec(Ndim))
      call Prec_calculate(E,Mprec,M,NCO,Ndim)

      allocate(Pk(Ndim)); Pk = 0.0D0
      call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

      allocate(Pmat(M,M),F2e(M,M),Fxc(M,M),Ftot(M,M))
      allocate(ApIA(Ndim),X(Ndim)); X = 0.0D0
      do iter=1,maxIter
         print*, "PCG_ITER:",iter

         ! CONVERT TRIAL VECTORS TO AO BASIS
         call VecToMat(Pk,Pmat,Coef,Ndim,NCO,M)

         ! CALCULATE TWO ELECTRON PART
         F2e = 0.0D0
         call g2g_calculate2e(Pmat,K4,cbas,1,F2e,1)
         F2e = 2.0D0 * F2e

         ! CALCULATE XC PART
         Fxc = 0.0D0
         call g2g_calculateg(Pmat,Fxc,2)

         ! OBTAIN FOCK TOTAL AND ADD TERM (Ea-Ei)Pk AND
         ! CHANGE BASIS AO -> MO
         call total_fock(F2e,Fxc,Ftot,M)
         call Ap_calculate(Ftot,Pk,Coef,E,ApIA,M,NCO,Nvirt,Ndim)
         
         ! CALCULATE ALPHA
         call Alpha_calc(Pk,ApIA,alpha,Ndim)

         ! NEW X VECTOR
         X = X + alpha * Pk

         ! NEW R VECTOR
         R = R - alpha * APIA

         ! CHECK CONVERGENCE
         call error(R,conv,Ndim)
         if (conv .eqv. .true.) exit
 
         ! GET NEW BETA AND Pk
         call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

      enddo ! ENDDO LOOP PCG

!     FORM RELAXED DENSITY OF EXCITED STATE
      call RelaxedDensity(X,Rho_urel,Coef,M,NCO,Ndim)
  
      deallocate(R,Mprec,Pk,Pmat,F2e,Fxc,Ftot,ApIA,X)
   end subroutine PCG_solve

   subroutine RelaxedDensity(Z,Rho_urel,C,M,NCO,N)
   use garcha_mod, only: RMM
      implicit none

      integer, intent(in) :: M, NCO, N
      real*8, intent(in) :: Z(N), C(M,M), Rho_urel(M,M)

      integer :: i, j, Nvirt, NCOc, pos, M2
      real*8, dimension(:,:), allocatable :: Rho_fund, Rho_exc, Rel_diff
      real*8, dimension(:,:), allocatable :: Zmo, Zao, RhoRMM

!     EXTRACT RHO FUND FROM RMM
      allocate(Rho_fund(M,M),RhoRMM(M,M))
      call spunpack_rho('L',M,RMM,Rho_fund)

!     CONVERT Z IN AO BASIS     
      allocate(Zmo(M,M)); Zmo = 0.0D0
      Nvirt = M - NCO
      NCOc = NCO + 1
      do i=1,NCO
      do j=1,Nvirt
         pos = (i - 1) * Nvirt + j
         Zmo(NCOc-i,NCO+j) = Z(pos)
      enddo
      enddo
   
      allocate(Zao(M,M))
      call matMOtomatAO(Zmo,Zao,C,M,1,.false.)
      deallocate(Zmo)

      allocate(Rel_diff(M,M))
      Rel_diff = Rho_urel + Zao

      allocate(Rho_exc(M,M))
      do i=1,M
      do j=1,M
         Rho_exc(i,j) = Rho_fund(i,j) + Rel_diff(i,j) + Rel_diff(j,i)
      enddo
      enddo

      print*, "Rho fundamental"
      do i=1,M
      do j=1,M
         print*, i,j,Rho_fund(i,j)
      enddo
      enddo

      print*, "relaxed density excited state"
      do i=1,M
      do j=1,M
         print*, i,j,Rho_exc(i,j)
      enddo
      enddo

!     SAVE EXCITED RHO INTO RMM
      M2 = M * 2
      do i=1,M
         RMM(i + (M2-i)*(i-1)/2) = Rho_exc(i,i)
         do j = i+1, M
            RMM(j + (M2-i)*(i-1)/2) = Rho_exc(i,j) * 2.0D0
         enddo
      enddo
   end subroutine RelaxedDensity

   subroutine error(V,convergence,N)
      implicit none

      integer, intent(in) :: N
      real*8, intent(in) :: V(N)
      logical, intent(inout) :: convergence
 
      integer :: i
      real*8 :: temp, tol
  
      tol = 1.0D-12
      temp = 0.0D0
      do i=1,N
        temp = temp + V(i)*V(i)
      enddo

      if( temp < tol ) convergence = .true.
   end subroutine error

   subroutine Alpha_calc(P,A,alpha,N)
      implicit none

      integer, intent(in) :: N
      real*8, intent(in) :: P(N), A(N)
      real*8, intent(out) :: alpha

      integer :: i
      real*8 :: temp
    
      temp = 0.0D0
      do i=1,N
        temp = temp + P(i) * A(i)
      enddo
      alpha = 1.0D0 / temp
   end subroutine Alpha_calc

   subroutine Ap_calculate(Fp,P,C,E,Ap,M,NCO,Nvirt,Ndim)
   use lr_data, only: Cocc_trans, Cvir

      implicit none

      integer, intent(in) :: M, NCO, Nvirt, Ndim
      real*8, intent(in) :: Fp(M,M), C(M,M), E(M), P(Ndim)
      real*8, intent(out) :: Ap(Ndim)

      integer :: i, j, NCOc, pos
      real*8, dimension(:,:), allocatable :: scratch

      allocate(scratch(M,Nvirt))
      scratch = matmul(Fp,Cvir)
      scratch = matmul(Cocc_trans,scratch)

!     FORM A*p IN MO BASIS
      NCOc = NCO + 1
      do i=1,NCO
      do j=1,Nvirt
         pos = (i-1) * Nvirt + j
         Ap(pos) = scratch(NCOc-i,j) + (E(NCO+j) - E(NCOc-i)) * P(pos)
      enddo
      enddo

      deallocate(scratch)
   end subroutine Ap_calculate

   subroutine VecToMat(Vec,Mat,Coef,Ndim,NCO,M)
   use lr_data, only: Coef_trans

      implicit none

      integer, intent(in) :: Ndim, NCO, M
      real*8, intent(in) :: Vec(Ndim), Coef(M,M)
      real*8, intent(out) :: Mat(M,M)

      integer :: row, col, NCOc, Nvirt, pos
     
      Nvirt = M - NCO
      NCOc = NCO + 1

      Mat = 0.0D0
      do row=1,NCO
      do col=1,Nvirt
         pos = (row-1) * Nvirt + col
         Mat(NCOc-row,NCO+col) = Vec(pos)
      enddo
      enddo
  
      Mat = matmul(Coef,Mat)
      Mat = matmul(Mat,Coef_trans)
   end subroutine VecToMat

   subroutine Pbeta_calc(R,M,beta,P,N)
      implicit none

      integer, intent(in) :: N
      real*8, intent(in) :: R(N), M(N)
      real*8, intent(inout) :: P(N)
      real*8, intent(out) :: beta

      integer :: i
      real*8 :: temp

      temp = 0.0D0
      do i=1,N
         temp = temp + R(i) * R(i) * M(i)
      enddo
      beta = 1.0D0 / temp

      do i=1,N
        P(i) = P(i) + beta * M(i) * R(i)
      enddo
   end subroutine Pbeta_calc

   subroutine Prec_calculate(Ener,Mprec,M,NCO,Ndim)
      implicit none
 
      integer, intent(in) :: M, NCO, Ndim
      real*8, intent(in) :: Ener(M)
      real*8, intent(out) :: Mprec(Ndim)

      integer :: i, j, Nvirt, NCOc, pos
      real*8 :: temp

      Nvirt = M - NCO
      NCOc = NCO + 1

      do i=1,NCO
      do j=1,Nvirt
         pos = (i-1) * Nvirt + j
         temp = Ener(j+NCO) - Ener(NCOc-i)
         temp = 1.0D0 / temp
         Mprec(pos) = temp
      enddo
      enddo
   end subroutine Prec_calculate

   subroutine total_fock(F1,F2,FT,M)
      implicit none

      integer, intent(in) :: M
      real*8, intent(in) :: F1(M,M)
      real*8, intent(inout) :: F2(M,M)
      real*8, intent(out) :: FT(M,M)

      integer :: i, j
      do i=1,M
      do j=i,M
         F2(j,i) = F2(i,j)
      enddo
      enddo
  
      FT = F1 + 2.0D0 * F2
   end subroutine total_fock

   subroutine basis_init(Coef,M,NCO,Nvirt)
!  This subroutine initializes the matrix needed for 
!  the change of basis in linear response and CPKS calculations.
      use lr_data, only: Coef_trans, Cocc, &
                        Cocc_trans, Cvir, Cvir_trans
      implicit none
      integer, intent(in) :: M, NCO, Nvirt
      real*8, intent(in) :: Coef(M,M)

      integer :: i, j

      allocate(Coef_trans(M,M),Cocc(M,NCO),Cocc_trans(NCO,M), &
               Cvir(M,Nvirt),Cvir_trans(Nvirt,M))

      do j=1,NCO
         Cocc(:,j) = Coef(:,j)
         Cocc_trans(j,:) = Coef(:,j)
      enddo
 
      do j=1,Nvirt
         Cvir(:,j) = Coef(:,NCO+j)
         Cvir_trans(j,:) = Coef(:,NCO+j)
      enddo

      do i=1,M
      do j=1,M
         Coef_trans(j,i) = Coef(i,j)
      enddo
      enddo
   end subroutine basis_init

   subroutine basis_deinit()
   use lr_data, only: Coef_trans, Cocc, &
                      Cocc_trans, Cvir, Cvir_trans
      implicit none
      deallocate(Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans)
   end subroutine basis_deinit
end module lrtddft
