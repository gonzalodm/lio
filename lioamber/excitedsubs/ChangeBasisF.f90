subroutine ChangeBasisF(FX,FT,Gxc,C,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
use lrdata, only: Cocc, Cocc_trans, Cvir, Cvir_trans

   implicit none

   integer, intent(in) :: M, Nvirt, NCO
   real*8, intent(in) :: C(M,M), Gxc(M,M), FX(M,M), FT(M,M)
   real*8, intent(out) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
   real*8, intent(out) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)

   integer :: i, j
   real*8, dimension(:,:), allocatable :: scratch

!  FORM FX IN BASIS VIRT X VIRT
   allocate(scratch(M,Nvirt))
   call multlr(FX,Cvir,scratch,M,M,Nvirt,1.0D0,0.0D0)
   call multlr(Cvir_trans,scratch,FXAB,Nvirt,M,Nvirt,1.0D0,0.0D0)
   deallocate(scratch)

!  FORM FX IN BASIS OCC X OCC
   allocate(scratch(M,NCO))
   call multlr(FX,Cocc,scratch,M,M,NCO,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,FXIJ,NCO,M,NCO,1.0D0,0.0D0)
   deallocate(scratch)

!  FORM FT IN BASIS OCC X VIR
   allocate(scratch(M,Nvirt))
   call multlr(FT,Cvir,scratch,M,M,Nvirt,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,FTIA,NCO,M,Nvirt,1.0D0,0.0D0)
   deallocate(scratch)

!  FORM GXC IN BASIS OCC X VIR
   allocate(scratch(M,Nvirt))
   call multlr(Gxc,Cvir,scratch,M,M,Nvirt,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,GXCIA,NCO,M,Nvirt,1.0D0,0.0D0)
   deallocate(scratch)
end subroutine ChangeBasisF
