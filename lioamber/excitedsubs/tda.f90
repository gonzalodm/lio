subroutine linear_response(MatCoef,VecEne,Xlr,deltaE,M,Nvirt,NCO,dim)
use lrdata, only: nstates, cbas, root, fitLR, Coef_trans, ppTDA
   implicit none

   integer, intent(in) :: M, Nvirt, NCO
   integer, intent(inout) :: dim
   real*8, intent(in) :: MatCoef(M,M), VecEne(M)
   real*8, intent(out) :: deltaE
   real*8, intent(out) :: Xlr(dim)

!  DAVIDSON OPTIONS
   integer :: maxIter, iter ! ITERATION
   integer :: vec_dim, iv, first_vec, newvec ! VECTOR INDEXES
   integer :: max_subs, Subdim ! DIMENSION SUBSPACE
   real*8, dimension(:,:), allocatable :: tvecMO,vmatAO,Fv,AX,H,eigvec
   real*8, dimension(:,:), allocatable :: RitzVec,ResMat
   real*8, dimension(:), allocatable :: eigval, val_old, Osc
   real*8, dimension(:,:,:), allocatable :: tmatMO,tmatAO,FM2,FXC
   integer :: calc_2elec
   logical :: conv

   integer :: i , j , k ! ELIMINAR

   calc_2elec = 0
   conv = .false.

   call g2g_timer_start('LINEAR RESPONSE')

   if (.not. ppTDA) then
      print*, ""
      print*,"======================================="
      print*,"         LINEAR RESPONSE - TDA"
      print*,"======================================="
      print*, ""
   else
      print*, ""
      print*,"======================================="
      print*,"       LINEAR RESPONSE - ppTDA"
      print*,"======================================="
      print*, ""
   endif

!  DIMENSION FOR ppTDA
   if ( ppTDA ) then
      dim = Nvirt * (Nvirt-1) / 2
   endif

!  CHECK INPUT VARIABLES FOR TDA
   if (nstates > dim) then
      print*, "NUMBER OF EXCITED STATES THAT YOU WANT IS BIGGER &
               & THAN DIMENSION MATRIX"
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
     allocate(tmatMO(M,M,vec_dim))
     call vecMOtomatMO(tvecMO,tmatMO,M,NCO,Nvirt, &
                       Subdim,vec_dim,first_vec,dim)

!    CHANGE BASIS MO -> AO FOR EACH VECTOR
     if(allocated(tmatAO)) deallocate(tmatAO)
        allocate(tmatAO(M,M,vec_dim))
     call matMOtomatAO(tmatMO,tmatAO,MatCoef,M,vec_dim,.true.)

!    CALCULATE 2E INTEGRALS
     allocate(FM2(M,M,vec_dim)); FM2 = 0.0D0
     call g2g_timer_start('2E integrals of vectrs')
     if (.not. ppTDA) then ! TDA
        if (.not. fitLR) then
           call g2g_calculate2e(tmatAO,cbas,vec_dim,FM2,calc_2elec)
           calc_2elec = 1 ! TURN OFF CALCULATE INTEGRALS
        else
           call calc2eFITT(tmatAO,FM2,vec_dim,M)
        endif
     else ! ppTDA
        call g2g_calcpptda(tmatAO,cbas,vec_dim,FM2)
       !do k=1,vec_dim
       !do i=1,M
       !do j=1,M
       !   print*, k,i,j,FM2(i,j,k)
       !enddo
       !enddo
       !enddo
       !stop
      
     endif
        
     call g2g_timer_stop('2E integrals of vectrs')

! THIS ONLY WORK IN TDA CALCULATION
     if (.not. ppTDA) then
!       CALCULATE DFT INTEGRALS
        allocate(FXC(M,M,vec_dim))
        call g2g_timer_start('Dft integrals of vectors')
        do iv=1,vec_dim ! LOOPS VECTORS INSIDE OF CYCLE
           call multlr(tmatMO(:,:,iv),Coef_trans,vmatAO,M,M,M,1.0D0,0.0D0)
           call g2g_calculatedft(vmatAO,MatCoef,Fv,NCO)
           FXC(:,:,iv) = Fv
           Fv = 0.0D0
        enddo ! END LOOPS VECTORS
        deallocate(tmatMO)
        call g2g_timer_stop('Dft integrals of vectors')

!       NOW FM2 CONTAIN THE 2e- AND DFT PARTS
        FM2 = FM2 + FXC 
        deallocate(FXC)
     endif

!    WE OBTAIN AX (NDIM x NVEC)  ! toy aqui
     call MtoIANV(FM2,MatCoef,AX,M,NCO,dim, &
                  Subdim,vec_dim,first_vec)
     deallocate(FM2)

!    ADD (Ea - Ei)*Xia TO AX
     call addInt(AX,VecEne,tvecMO,dim,M, &
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
       call OscStr(RitzVec,eigval,MatCoef,Osc,M,NCO,Nvirt,dim,nstates)
       call PrintResults(RitzVec,eigval,Osc,dim,nstates,M,NCO)
       exit
     endif

!    OBTAIN RESIDUAL VECTORS
     call residual(eigval,eigvec,RitzVec,AX,ResMat, &
                   dim,Subdim,nstates)

!    CHECK CONVERGENCE AND ADD NEW VECTORS
     conv = .false. ; newvec = 0
     call new_vectors(ResMat,eigval,VecEne,tvecMO,val_old, &
                      dim,Subdim,nstates,M,NCO,newvec,conv)

     if(conv .eqv. .true.) then
       write(*,*) ""
       write(*,"(1X,A,I2,1X,A)") "CONVERGED IN:",iter,"ITERATIONS"
       call OscStr(RitzVec,eigval,MatCoef,Osc,M,NCO,Nvirt,dim,nstates)
       call PrintResults(RitzVec,eigval,Osc,dim,nstates,M,NCO)
       exit
     else
!      ACTUALIZATION OF VECTORS INDEXES
       first_vec = Subdim 
       Subdim = Subdim + newvec
       vec_dim = newvec
     endif
   enddo !END DAVIDSON

   deallocate(vmatAO,Fv,ResMat,val_old,Osc,eigvec,&
              AX,tvecMO,tmatAO)

   call g2g_timer_stop('LINEAR RESPONSE')

   ! Save outputs
   if (root > 0 ) then
     Xlr = RitzVec(:,root)  
     deltaE = eigval(root)
   else
     Xlr = 0.0D0
     deltaE = 0.0D0
   endif

   deallocate(RitzVec,eigval)
end subroutine linear_response
