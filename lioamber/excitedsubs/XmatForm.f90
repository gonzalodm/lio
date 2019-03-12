subroutine XmatForm(Vec,Coef,Mat,Ndim,NCO,Nvirt,M)
use lrdata, only: Coef_trans
   implicit none

   integer, intent(in) :: Ndim, NCO, Nvirt, M
   real*8, intent(in) :: Vec(Ndim), Coef(M,M)
   real*8, intent(out) :: Mat(M,M)

   integer :: NCOc, row, col, pos
   real*8, dimension(:,:), allocatable ::SCR

   NCOc = NCO + 1
   do row=1,NCO
   do col=1,Nvirt
     pos = (row-1) * Nvirt + col
     Mat(NCOc-row,NCO+col) = Vec(pos)
   enddo
   enddo

   allocate(SCR(M,M))
   call multlr(Coef,Mat,SCR,M,M,M,1.0D0,0.0D0)
   call multlr(SCR,Coef_trans,Mat,M,M,M,1.0D0,0.0D0)
   deallocate(SCR)
end subroutine XmatForm
