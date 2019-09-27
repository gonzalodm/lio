subroutine open_OscStr(Xa,Xb,Ene,Ca,Cb,OsSt, &
                       M,NCOa,NCOb,Ndim,Nstat)
use lrdata, only: second_LR, doing_SLR, Tdip_flr, Tdip_slr, &
                  state_LR
   implicit none

   integer, intent(in) :: M, NCOa, NCOb, Ndim, Nstat
   real*8, intent(in) :: Xa(Ndim,Nstat), Xb(Ndim,Nstat), Ene(Nstat)
   real*8, intent(in) :: Ca(M,M), Cb(M,M)
   real*8, intent(out) :: OsSt(Nstat)

   integer :: Nvirta, Nvirtb
   real*8, dimension(:,:), allocatable :: Tdip
   real*8, dimension(:,:,:), allocatable :: TdensMOa,TdensMOb
   real*8, dimension(:,:,:), allocatable :: TdensAOa,TdensAOb

   Nvirta = M - NCOa
   Nvirtb = M - NCOb

   allocate(Tdip(Nstat,3),TdensMOa(M,M,Nstat),TdensMOb(M,M,Nstat))
   allocate(TdensAOa(M,M,Nstat),TdensAOb(M,M,Nstat))

!  FORM TRANSITION DENSITY IN MO BASIS
   ! ALPHA
   call vecMOtomatMO(Xa,TdensMOa,M,NCOa,Nvirta,Nstat,Nstat,0,Ndim)
   ! BETA
   call vecMOtomatMO(Xb,TdensMOb,M,NCOb,Nvirtb,Nstat,Nstat,0,Ndim)

!  CHANGE BASIS MO -> AO: ALPHA AND BETA
   call open_matMOtomatAO(tdensMOa,tdensMOb,tdensAOa,tdensAOb,&
        Ca,Cb,M,Nstat,.true.)
   deallocate(tdensMOa,tdensMOb)

!  CALCULATE TRANSITION DIPOLE: ANPHA AND BETA
   TdensAOa = TdensAOa + TdensAOb
   call TransDipole(TdensAOa,Tdip,M,Nstat)
   !Tdip = Tdip * 2.0D0 * dsqrt(2.0D0)
   deallocate(TdensAOa,TdensAOb)

!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat)
   OsSt = OsSt * 0.50D0
   
!  SAVE TRANSITION DIPOLE TO SECOND LINEAR RESPONSE
   if (second_LR) then
      print*, ""
      print*, "Saving unperturbed transition dipole"
      print*, ""
      if(allocated(Tdip_flr)) deallocate(Tdip_flr)
      allocate(Tdip_flr(Nstat,3))
      Tdip_flr = Tdip
   endif
   if (doing_SLR) then
      print*, ""
      print*, "Saving perturbed transition dipole"
      print*, ""
      if(allocated(Tdip_slr)) deallocate(Tdip_slr)
      allocate(Tdip_slr(Nstat,3))
      Tdip_slr = Tdip
   endif

   deallocate(Tdip)
end subroutine open_OscStr
