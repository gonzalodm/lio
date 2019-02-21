subroutine Pbeta_calc(R,M,beta,P,N)
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: R(N), M(N)
   real*8, intent(inout) :: P(N)
   real*8, intent(out) :: beta

   integer :: i
   real*8 :: temp

   temp = 0.0D0
   do i=1,N
      temp = temp + R(i) * R(i) * M(i)
   enddo
   beta = 1.0D0 / temp

   do i=1,N
     P(i) = P(i) + beta * M(i) * R(i)
   enddo
end subroutine Pbeta_calc
