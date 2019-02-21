subroutine total_fock(F1,F2,FT,M)
   implicit none

   integer, intent(in) :: M
   real*8, intent(in) :: F1(M,M)
   real*8, intent(inout) :: F2(M,M)
   real*8, intent(out) :: FT(M,M)

   integer :: i, j
   do i=1,M
   do j=i,M
      F2(j,i) = F2(i,j)
   enddo
   enddo

   FT = F1 + 2.0D0 * F2
end subroutine total_fock
