subroutine ObtainOsc(dip,E,O,N)
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: dip(N,3), E(N)
   real*8, intent(out) :: O(N)

   integer :: i, j
   real*8 :: dostres = 2.0D0 / 3.0D0
   real*8 :: temp = 0.0D0

   do i=1,N
   do j=1,3
      temp = temp + dip(i,j) * dip(i,j) * E(i) * dostres
   enddo
   O(i) = temp; temp = 0.0D0
   enddo
end subroutine ObtainOsc
