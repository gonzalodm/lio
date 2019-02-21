subroutine VecToMat(Vec,Mat,Coef,Ndim,NCO,M)
use lrdata, only: Coef_trans

   implicit none

   integer, intent(in) :: Ndim, NCO, M
   real*8, intent(in) :: Vec(Ndim), Coef(M,M)
   real*8, intent(out) :: Mat(M,M)

   integer :: row, col, NCOc, Nvirt, pos
   real*8, dimension(:,:), allocatable :: scr

   Nvirt = M - NCO
   NCOc = NCO + 1

   Mat = 0.0D0
   do row=1,NCO
   do col=1,Nvirt
      pos = (row-1) * Nvirt + col
      Mat(NCOc-row,NCO+col) = Vec(pos)
   enddo
   enddo

   allocate(scr(M,M))
   call multlr(Coef,Mat,scr,M,M,M,1.0D0,0.0D0)
   call multlr(scr,Coef_trans,Mat,M,M,M,1.0D0,0.0D0)
   deallocate(scr)
end subroutine VecToMat
