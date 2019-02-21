subroutine formred(TmatMO,C,VecMOAO,M,dim,numvec,iv)
use lrdata, only: Coef_trans

  implicit none

  integer, intent(in) :: M, dim, iv, numvec
  real*8, intent(in) :: TmatMO(M,M,numvec)
  real*8, intent(in) :: C(M,M)
  real*8, intent(out) :: VecMOAO(M,M)

  integer :: i, j

  call multlr(TmatMO(:,:,iv),Coef_trans,vecMOAO,M,M,M,1.0D0,0.0D0)
end subroutine formred
