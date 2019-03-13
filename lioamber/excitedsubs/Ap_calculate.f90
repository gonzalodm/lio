subroutine Ap_calculate(Fp,P,E,Ap,M,NCO,Nvirt,Ndim)
use lrdata, only: Cocc_trans, Cvir

   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Ndim
   real*8, intent(in) :: Fp(M,M), E(M), P(Ndim)
   real*8, intent(out) :: Ap(Ndim)

   integer :: i, j, NCOc, pos
   real*8, dimension(:,:), allocatable :: scrA, scrB

   allocate(scrA(M,Nvirt),scrB(NCO,Nvirt))
   call multlr(Fp,Cvir,scrA,M,M,Nvirt,1.0D0,0.0D0)
   call multlr(Cocc_trans,scrA,scrB,NCO,M,Nvirt,1.0D0,0.0D0)

!  FORM A*p IN MO BASIS
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,Nvirt
      pos = (i-1) * Nvirt + j
      Ap(pos) = scrB(NCOc-i,j) + ( E(NCO+j) - E(NCOc-i) )*P(pos)
   enddo
   enddo

   deallocate(scrA,scrB)
end subroutine Ap_calculate
