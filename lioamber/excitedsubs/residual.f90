subroutine residual(Val,Vec,R,A,W,Ndim,Sdim,Nstat)
   implicit none

   integer, intent(in) :: Ndim, Sdim, Nstat
   real*8, intent(in) :: Val(Nstat), Vec(Sdim,Nstat), R(Ndim,Nstat)
   real*8, intent(in) :: A(Ndim,Sdim)
   real*8, intent(out) :: W(Ndim,Nstat)

   integer :: i
   real*8, dimension(:), allocatable :: temp

   allocate(temp(Ndim))
   do i=1,Nstat
     call multlr(A,Vec(:,i),temp,Ndim,Sdim,1,1.0D0,0.0D0)
     W(1:Ndim,i) = temp - Val(i) * R(1:Ndim,i)
   enddo
   deallocate(temp)
end subroutine residual
