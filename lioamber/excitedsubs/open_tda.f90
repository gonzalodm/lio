subroutine open_linear_response(Ca,Cb,Ea,Eb,M,NCOa,NCOb)
use lrdata, only: nstates
   implicit none

   integer, intent(in) :: M, NCOa, NCOb
   real*8, intent(in) :: Ca(M,M), Cb(M,M)
   real*8, intent(in) :: Ea(M), Eb(M)

   integer :: i, j, k
   integer :: iter, calc_2elec
   integer :: Ndima, Ndimb, Ndimt, Ndim, Nvirta, Nvirtb
   integer :: vec_dim, newvec, first_vec, Subdim, maxIter, max_subs
   real*8, dimension(:), allocatable :: val_old, Osc
   real*8, dimension(:,:), allocatable :: AX, tvecMOa, tvecMOb
   real*8, dimension(:,:,:), allocatable :: tmatMOa, tmatMOb, tmatAOa, tmatAOb
   real*8, dimension(:,:,:), allocatable :: FM2a, FM2b

   calc_2elec = 0
   ! Set occ and virtual and dimension of full matrix
   Nvirta = M - NCOa
   Nvirtb = M - NCOb
   Ndima = NCOa*Nvirta
   Ndimb = NCOb*Nvirtb
   Ndimt = Ndima + Ndimb
   Ndim = max(Ndima,Ndimb)
   print*, "M,Ndim,Ndimt,NCOa,NCOb",M,Ndim,Ndimt,NCOa,NCOb

   print*, ""
   print*,"======================================="
   print*,"      OPEN LINEAR RESPONSE - TDA"
   print*,"======================================="
   print*, ""

   !  CHECK INPUT VARIABLES
   if (nstates > Ndim) then
      print*, "NUMBER OF EXCITED STATES THAT YOU WANT IS BIGGER &
               & THAN DIMENSION MATRIX"
      nstates = Ndim
      print*, "WE CALCULATED", nstates
   endif

   ! DAVIDSON INITIALIZATION
   allocate(val_old(nstates),Osc(nstates))

   val_old = 1.0D0
   vec_dim = 4 * nstates
   first_vec = 0
   newvec = 0

   if(vec_dim >= Ndim) then
      vec_dim = Ndim
      Subdim = Ndim
      maxIter = 1
      max_subs = Ndim
      allocate(AX(Ndim,Ndim),tvecMOa(Ndim,Ndim),tvecMOb(Ndim,Ndim))
   else
      max_subs = 200
      if (max_subs > Ndim) max_subs = Ndim
      maxIter = 50
      Subdim = vec_dim
      allocate(AX(Ndim,max_subs),tvecMOa(Ndim,max_subs),&
               tvecMOb(Ndim,max_subs))
   endif

   ! INIT TRIAL VECTORS
   call vec_init(tvecMOa,Ndim,vec_dim)
   call vec_init(tvecMOb,Ndim,vec_dim)

!  PRINT INFORMATION
   write(*,"(1X,A,10X,I2,2X,I2,2X,I2,2X,I2)") "NCOA, NCOB, NVIRTA, NVIRTB",NCOa,NCOb,Nvirta,Nvirtb
   write(*,"(1X,A,10X,I2,2X,I2,2X,I2)") "NDIMA, NDIMB, M", Ndima, Ndimb, M
   write(*,"(1X,A,11X,I3)") "DIMENSION OF FULL MATRIX",Ndim
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
     ! alpha
     allocate(tmatMOa(M,M,vec_dim))
     call vecMOtomatMO(tvecMOa,tmatMOa,M,NCOa,Nvirta, &
                       Subdim,vec_dim,first_vec,Ndim)
     ! beta
     allocate(tmatMOb(M,M,vec_dim))
     call vecMOtomatMO(tvecMOb,tmatMOb,M,NCOb,Nvirtb, &
                       Subdim,vec_dim,first_vec,Ndim)

!    print*, "MatMO delos vec, alpha, beta"
!     do k=1,vec_dim
!     print*, "vector", k
!    do i=1,M
!    do j=1,M
!       print*, i,j,tmatMOa(i,j,1),tmatMOb(i,j,1)
!    enddo
!    enddo
!     enddo

!    CHANGE BASIS MO -> AO FOR EACH VECTOR
     ! alpha
     if(allocated(tmatAOa)) deallocate(tmatAOa)
        allocate(tmatAOa(M,M,vec_dim))
     ! beta
     if(allocated(tmatAOb)) deallocate(tmatAOb)
        allocate(tmatAOb(M,M,vec_dim))
     call open_matMOtomatAO(tmatMOa,tmatMOb,tmatAOa,tmatAOb,&
          Ca,Cb,M,vec_dim,.true.)
!    print*, "MatAO delos vec, alpha, beta"
!     do k=1,vec_dim
!     print*, "vector", k
!    do i=1,M
!    do j=1,M
!       print*, i,j,tmatAOa(i,j,1),tmatAOb(i,j,1)
!    enddo
!    enddo
!     enddo

!    CALCULATE 2E INTEGRALS
     allocate(FM2a(M,M,vec_dim)); FM2a = 0.0D0
     allocate(FM2b(M,M,vec_dim)); FM2b = 0.0D0
     call g2g_timer_start('2E integrals of vectrs')
     if (.not. fitLR) then
        call g2g_open_calculate2e(tmatAOa,tmatAOb,cbas,vec_dim,
                                  FM2a,FM2b,calc_2elec)
        calc_2elec = 1 ! TURN OFF CALCULATE INTEGRALS
     else
        print*, "NOT FITTIN IN OPEN SHELL LINEAR RESPONSE"
        stop
        !call calc2eFITT(tmatAO,FM2,vec_dim,M)
     endif
     call g2g_timer_stop('2E integrals of vectrs')





     

   enddo ! END DAVIDSON CYCLE
     
end subroutine open_linear_response
