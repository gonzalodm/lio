subroutine open_linear_response(Ca,Cb,Ea,Eb,M,NCOa,NCOb)
use lrdata, only: nstates, cbas, fitLR, Coef_trans, Coef_transB
   implicit none

   integer, intent(in) :: M, NCOa, NCOb
   real*8, intent(in) :: Ca(M,M), Cb(M,M)
   real*8, intent(in) :: Ea(M), Eb(M)

   integer :: i, j, k, iv
   integer :: iter, calc_2elec
   integer :: Ndima, Ndimb, Ndimt, Ndim, Nvirta, Nvirtb
   integer :: vec_dim, newvec, first_vec, Subdim, maxIter, max_subs
   real*8, dimension(:), allocatable :: val_old, Osc, eigval
   real*8, dimension(:,:), allocatable :: AX_a, AX_b, tvecMOa, tvecMOb, vmatAOa, vmatAOb
   real*8, dimension(:,:), allocatable :: Fva, Fvb, H, eigvec, RitzVecA, RitzVecB
   real*8, dimension(:,:), allocatable :: ResMatA, ResMatB
   real*8, dimension(:,:,:), allocatable :: tmatMOa, tmatMOb, tmatAOa, tmatAOb
   real*8, dimension(:,:,:), allocatable :: FM2a, FM2b, FXC_a, FXC_b
   logical :: conv

   calc_2elec = 0

   ! SET DIMENSION OF ALPHA AND BETA MATRIX AND DIMENSION OF FULL MATRIX
   Nvirta = M - NCOa
   Nvirtb = M - NCOb
   Ndima = NCOa * Nvirta
   Ndimb = NCOb * Nvirtb
   Ndimt = Ndima + Ndimb
   Ndim  = max(Ndima,Ndimb)

   print*, ""
   print*,"======================================="
   print*,"      OPEN LINEAR RESPONSE - TDA"
   print*,"======================================="
   print*, ""

   !  CHECK INPUT VARIABLES
   if (nstates > Ndimt) then
      print*, "NUMBER OF EXCITED STATES THAT YOU WANT IS BIGGER &
               & THAN DIMENSION MATRIX"
      nstates = Ndimt
      print*, "WE CALCULATED", nstates
   endif

   ! DAVIDSON INITIALIZATION
   allocate(val_old(nstates),Osc(nstates))

   val_old = 1.0D0
   vec_dim = 4 * nstates
   first_vec = 0
   newvec = 0

   if(vec_dim >= Ndimt) then
      vec_dim = Ndimt
      Subdim = Ndimt
      maxIter = 1
      max_subs = Ndimt
      allocate(tvecMOa(Ndim,Ndimt),tvecMOb(Ndim,Ndimt))
      allocate(AX_a(Ndim,Ndimt),AX_b(Ndim,Ndimt))
   else
      max_subs = 200
      if (max_subs > Ndimt) max_subs = Ndimt
      maxIter = 50
      Subdim = vec_dim
      allocate(tvecMOa(Ndim,max_subs),tvecMOb(Ndim,max_subs))
      allocate(AX_a(Ndim,max_subs),AX_b(Ndim,max_subs))
   endif

   allocate(vmatAOa(M,M),vmatAOb(M,M))
   allocate(Fva(M,M),Fvb(M,M))
   allocate(RitzVecA(Ndim,nstates),RitzVecB(Ndim,nstates))
   allocate(ResMatA(Ndim,nstates),ResMatB(Ndim,nstates))
   Fva = 0.0D0
   Fvb = 0.0D0
   RitzVecA = 0.0D0; RitzVecB = 0.0D0

   ! INIT TRIAL VECTORS
   call open_vec_init(tvecMOa,tvecMOb,Ea,Eb,M,NCOa,NCOb,Ndim,vec_dim)

!  PRINT INFORMATION
   write(*,"(1X,A,10X,I2,2X,I2,2X,I2,2X,I2)") "NCOA, NCOB, NVIRTA, NVIRTB",NCOa,NCOb,Nvirta,Nvirtb
   write(*,"(1X,A,10X,I2,2X,I2,2X,I2)") "NDIMA, NDIMB, M", Ndima, Ndimb, M
   write(*,"(1X,A,11X,I3)") "DIMENSION OF FULL MATRIX",Ndimt
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

!    CHANGE BASIS MO -> AO FOR EACH VECTOR
     ! alpha
     if(allocated(tmatAOa)) deallocate(tmatAOa)
        allocate(tmatAOa(M,M,vec_dim))
     ! beta
     if(allocated(tmatAOb)) deallocate(tmatAOb)
        allocate(tmatAOb(M,M,vec_dim))
     call open_matMOtomatAO(tmatMOa,tmatMOb,tmatAOa,tmatAOb,&
          Ca,Cb,M,vec_dim,.true.)

!    CALCULATE 2E INTEGRALS
     allocate(FM2a(M,M,vec_dim)); FM2a = 0.0D0
     allocate(FM2b(M,M,vec_dim)); FM2b = 0.0D0
     call g2g_timer_start('2E integrals of vectrs')
     if (.not. fitLR) then
        call g2g_open_calculate2e(tmatAOa,tmatAOb,cbas,vec_dim, &
                                  FM2a,FM2b,calc_2elec)
        calc_2elec = 1 ! TURN OFF CALCULATE INTEGRALS
     else
        call open_calc2eFITT(tmatAOa,tmatAOb,FM2a,FM2b,vec_dim,M)
     endif
     call g2g_timer_stop('2E integrals of vectrs')
    
!    CALCULATE DFT INTEGRALS
     call g2g_timer_start('Dft integrals of vectors')
     allocate(FXC_a(M,M,vec_dim),FXC_b(M,M,vec_dim))
     do iv=1,vec_dim ! LOOPS VECTORS INSIDE OF ITERATION
        ! beta
        call multlr(tmatMOb(:,:,iv),Coef_transB,vmatAOb,M,M,M,1.0D0,0.0D0)
        ! alpha
        call multlr(tmatMOa(:,:,iv),Coef_trans,vmatAOa,M,M,M,1.0D0,0.0D0)
        call g2g_open_calculatedft(vmatAOa,vmatAOb,Ca,Cb,Fva,Fvb,NCOa,NCOb)
        FXC_a(:,:,iv) = Fva
        FXC_b(:,:,iv) = Fvb
        Fva = 0.0D0
        Fvb = 0.0D0
     enddo ! END LOOPS VECTORS
     deallocate(tmatMOa,tmatMOb)
     call g2g_timer_stop('Dft integrals of vectors')

!    OBTAIN TOTAL FOCK: ALPHA AND BETA
     FM2a = FM2a + FXC_a
     FM2b = FM2b + FXC_b
     deallocate(FXC_a,FXC_b)
     call open_MtoIANV(FM2a,FM2b,Ca,Cb,AX_a,AX_b,M,NCOa,NCOb,&
                       Ndim,Subdim,vec_dim,first_vec)
     deallocate(FM2a,FM2b)

!    ADD ENERGY
     ! alpha 
     call addInt(AX_a,Ea,tvecMOa,Ndim,M,Subdim,NCOa,vec_dim,first_vec)
     ! beta
     call addInt(AX_b,Eb,tvecMOb,Ndim,M,Subdim,NCOb,vec_dim,first_vec)

!    WE OBTAIN SUBSPACE MATRIX
     allocate(H(Subdim,Subdim))
     call open_subspaceMat(AX_a,AX_b,tvecMOa,tvecMOb,H,&
                           Ndima,Ndimb,Ndim,Subdim)

!    DIAGONALIZATION OF SUBSPACE MATRIX
     if(allocated(eigvec)) deallocate(eigvec)
     if(allocated(eigval)) deallocate(eigval)
       allocate(eigvec(Subdim,Subdim),eigval(Subdim))
     call diagonH(H,Subdim,eigval,eigvec)
     deallocate(H)

!    CHANGE SUBSPACE TO FULL SPACE - OBTAIN RITZ VECTOR
     call open_RitzObtain(tvecMOa,tvecMOb,eigvec,RitzVecA,RitzVecB, &
                          Ndim,Ndima,Ndimb,Subdim,nstates)
     
!    OBTAIN RESIDUAL VECTORS
     call open_residual(eigval,eigvec,RitzVecA,RitzVecB,AX_a,AX_b, &
                        ResMatA,ResMatB,Ndima,Ndimb,Ndim,Subdim,nstates)

!    CHECK CONVERGENCE AND ADD NEW VECTORS
     conv = .false.; newvec = 0
     call open_new_vectors(ResMatA,ResMatb,eigval,Ea,Eb,tvecMOa,tvecMOb, &
                           val_old,Ndima,Ndimb,Ndim,Subdim,nstates,M,NCOa, &
                           NCOb,newvec,conv)

     if(conv .eqv. .true.) then
       write(*,*) ""
       write(*,"(1X,A,I2,1X,A)") "CONVERGED IN:",iter,"ITERATIONS"
!      OSCILLATOR STRENGTH CALCULATION
       call open_OscStr(RitzVecA,RitzVecB,eigval,Ca,Cb,Osc, &
                       M,NCOa,NCOb,Ndim,nstates)

!      PRINT RESULTS
       call open_PrintResults(RitzVecA,RitzVecB,eigval,Osc, &
                              Ndim,nstates,M,NCOa,NCOb)
       exit
     else
!      ACTUALIZATION OF VECTORS INDEXES
       first_vec = Subdim
       Subdim = Subdim + newvec
       vec_dim = newvec
     endif
   enddo ! END DAVIDSON CYCLE

   deallocate(vmatAOa,vmatAOb,Fva,Fvb,ResMatA,ResMatB,val_old,Osc,eigvec,&
              AX_a,AX_b,tvecMOa,tvecMOb,tmatAOa,tmatAOb)

end subroutine open_linear_response
