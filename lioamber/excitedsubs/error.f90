subroutine error(V,convergence,N,iter)
  implicit none

  integer, intent(in) :: N, iter
  real*8, intent(in) :: V(N)
  logical, intent(inout) :: convergence

  integer :: i
  real*8 :: temp, tol

  tol = 1.0D-12
  temp = 0.0D0
  do i=1,N
    temp = temp + V(i)*V(i)
  enddo

  write(*,8070) iter, temp, tol
  print*, ""
  if( temp < tol ) convergence = .true.

8070  FORMAT(1X,"Iteration = ", I2,1X,&
      "Error (crit) =",ES9.2," (",ES9.2,")")
end subroutine error
