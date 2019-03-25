subroutine MtoIANV(F,C,A,M,NCO,Ndim,Sdim,Nvec,V1)
use lrdata, only: Cocc

   implicit none

   integer, intent(in) :: M, Ndim, Nvec, NCO, V1, Sdim
   real*8, intent(in) :: F(M,M,Nvec), C(M,M)
   real*8, intent(inout) :: A(Ndim,Sdim)

   integer :: i, j, k, iv, row, Nvirt, NCOc
   real*8, dimension(:,:), allocatable :: B

   real*8 :: temp

   Nvirt = M - NCO
   NCOc = NCO + 1
   allocate(B(NCO,M))

   temp = 0.0D0
   do iv=1,Nvec
     call multlr(F(:,:,iv),Cocc,B,M,M,NCO,1.0D0,0.0D0)

     do i=1,NCO
     do j=NCOc,M
       do k=1,M
         temp = temp + B(k,NCOc-i) * C(k,j)
       enddo
       row = (i-1) * Nvirt + (j-NCO)
       A(row,V1+iv) = temp
       temp = 0.0D0
     enddo
     enddo
   enddo

   deallocate(B)
end subroutine MtoIANV
