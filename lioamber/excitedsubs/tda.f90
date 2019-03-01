subroutine linear_response(MatCoef,VecEne)
use lrdata, only: Nvirt,dim,nstates,eigenval,&
                   eigenvec,cbas,root,FCA,nfo,nfv,&
                   NCOlr, Mlr, fitLR, EneSCF
use garcha_mod, only: NCO, M, c, a

   implicit none

   real*8, intent(in) :: MatCoef(M,M)
   real*8, intent(in) :: VecEne(M)

!  COUNTERS
   integer :: i, j
!  IF FCA IS USED
   real*8, dimension(:), allocatable :: Ene_LR
   real*8, dimension(:,:), allocatable :: Coef_LR
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
      allocate(Coef_LR(M,M),Ene_LR(M),EneSCF(M))
      Coef_LR = MatCoef
      Ene_LR = VecEne
      EneSCF = VecEne ! para forces
   endif

!  INITIALIZATION OF MATRIX NEEDED FOR CHANGE BASIS
   call basis_initLR(Coef_LR,M,NCO,Nvirt)

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

!  For forces
   allocate(eigenval(nstates))

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
   allocate(vmatAO(M,M),Fv(M,M))
   allocate(RitzVec(dim,nstates),ResMat(dim,nstates))
   vmatAO = 0.0D0; Fv = 0.0D0; RitzVec = 0.0D0

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
     call g2g_timer_start('2E integrals of vectrs')
     if (.not. fitLR) then
        call g2g_calculate2e(tmatAO,cbas,vec_dim,FM2,calc_2elec)
        calc_2elec = 1 ! TURN OFF CALCULATE INTEGRALS
     else
        call calc2eFITT(tmatAO,FM2,vec_dim,M)
     endif
     call g2g_timer_stop('2E integrals of vectrs')

!    CALCULATE DFT INTEGRALS
     allocate(FXC(M,M,vec_dim))
     call g2g_timer_start('Dft integrals of vectors')
     do iv=1,vec_dim ! LOOPS VECTORS INSIDE OF CYCLE
        call formred(tmatMO,Coef_LR,vmatAO,M,dim,vec_dim,iv) ! FORM T * C**T
        call g2g_calculatedft(vmatAO,Coef_LR,Fv,NCO)
        FXC(:,:,iv) = Fv
        Fv = 0.0D0
     enddo ! END LOOPS VECTORS
     call g2g_timer_stop('Dft integrals of vectors')

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
       write(*,*) ""
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

!  Copy energies for forces
   do i=1,nstates
      eigenval(i) = eigval(i)
   enddo

   deallocate(vmatAO,Fv,ResMat,val_old,Osc,eigvec,&
              eigval,AX,tvecMO,tmatAO)

   call g2g_timer_stop('LINEAR RESPONSE')

   if (root > 0 ) then
     call Zvector(Coef_LR,Ene_LR,RitzVec(:,root),NCO,M,dim)
   endif
   
   deallocate(RitzVec,Coef_LR,Ene_LR)
   call basis_deinitLR()
end subroutine linear_response
