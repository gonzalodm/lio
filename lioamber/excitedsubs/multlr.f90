subroutine multlr(A,B,C,M,K,N,alpha,beta)
   implicit none

   integer, intent(in) :: M, K, N
   real*8, intent(in) :: A(M,K), B(K,N)
   real*8, intent(in) :: alpha, beta
   real*8, intent(inout) :: C(M,N)

   call dgemm('N','N',M,N,K,alpha,A,M,B,K,beta,C,M)
end subroutine multlr
