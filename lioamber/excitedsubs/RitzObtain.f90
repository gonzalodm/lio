subroutine RitzObtain(V,Vsmall,R,Ndim,Sdim,Nstat)
   implicit none

   integer, intent(in) :: Ndim, Sdim, Nstat
   real*8, intent(in) :: V(Ndim,Sdim), Vsmall(Sdim,Nstat)
   real*8, intent(out) :: R(Ndim,Nstat)

   call multlr(V,Vsmall,R,Ndim,Sdim,Nstat,1.0D0,0.0D0)
end subroutine RitzObtain
