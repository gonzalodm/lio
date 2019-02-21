subroutine norma(V,N,norm)
  implicit none

  integer, intent(in) :: N
  real*8, intent(in) :: V(N)
  real*8,intent(out) :: norm

  integer :: i

  norm = 0.0d0
  do i=1,N
    norm = norm + V(i)*V(i)
  enddo
end subroutine norma
