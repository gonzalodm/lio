subroutine subspaceMat(A,V,H,Ndim,Sdim)
   implicit none

   integer, intent(in) :: Ndim, Sdim
   real*8, intent(in) :: A(Ndim,Sdim), V(Ndim,Sdim)
   real*8, intent(out) :: H(Sdim,Sdim)

   integer :: i, j
   real*8, dimension(:,:), allocatable :: VT

   allocate(VT(Sdim,Ndim))
   do i=1,Ndim
   do j=1,Sdim
      VT(j,i) = V(i,j)
   enddo
   enddo

   call multlr(VT,A,H,Sdim,Ndim,Sdim,1.0D0,0.0D0)
   deallocate(VT)
end subroutine subspaceMat
