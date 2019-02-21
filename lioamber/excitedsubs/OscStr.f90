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

!  FORM TRANSITION DENSITY IN MO BASIS
   call vecMOtomatMO(X,TdensMO,M,NCO,Nvirt,Nstat,Nstat,0,Ndim)
!  CHANGE BASIS MO -> AO
   call matMOtomatAO(TdensMO,TdensAO,Coef,M,Nstat,.true.)
   deallocate(TdensMO)
!  CALCULATE TRANSITION DIPOLE
   call TransDipole(TdensAO,Tdip,M,Nstat)
   deallocate(TdensAO)
!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat)
   deallocate(Tdip)
end subroutine OscStr
