subroutine diagonH(H,N,Val,Vec)
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: H(N,N)
   real*8, intent(out) :: Val(N),Vec(N,N)

   real*8, dimension(:), allocatable :: WORK
   integer, dimension(:), allocatable :: IWORK
   real*8 :: WI(N)
   integer :: i, j, info, LWORK, LIWORK

   Vec = H
   LWORK = -1
   LIWORK = -1

   allocate(WORK(1),IWORK(1))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   LWORK=WORK(1)
   LIWORK=IWORK(1)
   deallocate(WORK,IWORK); allocate(WORK(LWORK),IWORK(LIWORK))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   deallocate(WORK,IWORK)
end subroutine diagonH
