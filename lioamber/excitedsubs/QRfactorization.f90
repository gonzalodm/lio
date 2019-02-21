subroutine QRfactorization(A,N,M)
   implicit none

   integer, intent(in) :: N, M
   real*8, intent(inout) :: A(N,M)

   real*8,dimension(:),allocatable :: WORK, TAU
   integer :: i,j,LWORK,INFO,K
   real*8 :: norma, ortog

   allocate(WORK(1),TAU(M))
   call dgeqrf(N,M,A,N,TAU,WORK,-1,INFO)
   LWORK = WORK(1)
   deallocate(WORK); allocate(WORK(LWORK))
   call dgeqrf(N,M,A,N,TAU,WORK,LWORK,INFO)
   call dorgqr(N,M,M,A,N,TAU,WORK,LWORK,INFO)
end subroutine QRfactorization
