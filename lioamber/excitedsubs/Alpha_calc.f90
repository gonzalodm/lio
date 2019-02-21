subroutine Alpha_calc(P,A,alpha,N)
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: P(N), A(N)
   real*8, intent(out) :: alpha

   integer :: i
   real*8 :: temp

   temp = 0.0D0
   do i=1,N
     temp = temp + P(i) * A(i)
   enddo
   alpha = 1.0D0 / temp
end subroutine Alpha_calc
