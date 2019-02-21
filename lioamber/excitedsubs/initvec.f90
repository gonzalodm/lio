subroutine vec_init(Vec,N,vecnum)
   implicit none

   integer, intent(in) :: N, vecnum
   real*8, intent(out) :: Vec(N,vecnum)

   integer :: i

   Vec = 0.0D0
   do i=1,vecnum
      Vec(i,i) = 1.0D0
   enddo
end subroutine vec_init

