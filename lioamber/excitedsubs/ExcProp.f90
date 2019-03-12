subroutine ExcProp(C_scf,E_scf)
use lrdata, only: root, forEXC, nstates
use garcha_mod, only: M, NCO, writeforces, natom, doing_ehrenfest
use td_data, only: timedep
use maskrmm, only: rmmput_dens
   implicit none

   real*8, allocatable, intent(inout) :: C_scf(:,:), E_scf(:)

   integer :: NCOlr, Mlr, Nvirt, Ndim
   real*8 :: deltaE
   real*8, dimension(:), allocatable :: Xlr, Zvec, Qvec
   real*8, dimension(:,:), allocatable :: Punr, Trans, Gxc
   real*8, dimension(:,:), allocatable :: rhoEXC, Pdif

   ! for the moment fca is broken
   call fcaApp(C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)

   ! We form matrix for change basis
   call basis_initLR(C_scf,M,NCO,Nvirt)

   ! Linear Response Calculation
   allocate(Xlr(Ndim)) ! Vectors of LR
   Xlr = 0.0D0
   call linear_response(C_scf,E_scf,Xlr,deltaE,M,Nvirt,NCO,Ndim)
   Xlr = Xlr / dsqrt(2.0D0) ! We normalize to 1/2

   ! Z-Vector Method
   if ( root > 0 ) then
      if ( root > nstates ) then
         print*, "The excited state that you want was not &
         calculated in Linear Response. Please root <= nstates"
         stop
      endif

      ! Punr = Unrelaxed Difference Density Matrix in AO basis
      ! Trans = Transition Density Matrix in AO basis
      allocate(Punr(M,M),Trans(M,M))
      Punr = 0.0D0
      Trans = 0.0D0
      call GenerateDensities(Xlr,C_scf,Punr,Trans,&
                             M,Ndim,NCO,Nvirt)

      allocate(Zvec(Ndim),Qvec(Ndim),Gxc(M,M))
      Zvec = 0.0D0
      Qvec = 0.0D0
      Gxc = 0.0D0
      call Zvector(C_scf,E_scf,Xlr,Punr,Trans,&
                   Zvec,Qvec,Gxc,NCO,M,Ndim,Nvirt)

      ! OBTAIN DENSITIES
      allocate(rhoEXC(M,M),Pdif(M,M))
      rhoEXC = 0.0D0
      Pdif = 0.0D0
      ! rhoEXC = Total Excited Density Matrix
      ! Pdif = Relaxed Difference Density Matrix
      call RelaxedDensity(Zvec,Punr,C_scf,&
                          RhoEXC,Pdif,M,NCO,Ndim,Nvirt)
      deallocate(Punr)

      ! Put Excited state density in RMM
      if ( (timedep == 1 .or. doing_ehrenfest) &
            .and. (.not. writeforces) ) then
         print*, "--- Put Excited State Density Matrix in RMM ---"
         call rmmput_dens( rhoEXC )

      endif

      ! Forces in Excited State
      if ( writeforces ) then
         allocate(forEXC(natom,3))
         call forcesexc(rhoEXC,Pdif,Zvec,Trans,&
                        Qvec,Gxc,Xlr,C_scf,E_scf,deltaE,&
                        forEXC,M,Ndim,natom,NCO)
      endif 
      deallocate(Zvec,Qvec,Gxc,Trans,rhoEXC,Pdif)
   endif

   deallocate(Xlr)

   call basis_deinitLR()

end subroutine ExcProp
