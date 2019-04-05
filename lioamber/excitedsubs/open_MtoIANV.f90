subroutine open_MtoIANV(Fa,Fb,Ca,Cb,Aa,Ab,M,NCOa,NCOb,&
                       Ndim,Sdim,Nvec,V1)
use lrdata, only: Cocc_trans, Cocc_transB

   implicit none

   integer, intent(in) :: M, NCOa, NCOb, Ndim, Sdim, Nvec, V1
   real*8, intent(in) :: Fa(M,M,Nvec), Fb(M,M,Nvec), Ca(M,M), Cb(M,M)
   real*8, intent(inout) :: Aa(Ndim,Sdim), Ab(Ndim,Sdim)

   integer :: i, j, k, iv, row, Nvirt, NCOc
   real*8 :: temp
   real*8, dimension(:,:), allocatable :: B

   ! ALPHA
   Nvirt = M - NCOa
   NCOc = NCOa + 1
   allocate(B(NCOa,M))
   temp = 0.0D0
   do iv=1,Nvec
     call multlr(Cocc_trans,Fa(:,:,iv),B,NCOa,M,M,1.0D0,0.0D0)
     do i=1,NCOa
     do j=NCOc,M
       do k=1,M
         temp = temp + B(NCOc-i,k) * Ca(k,j)
       enddo
       row = (i-1) * Nvirt + (j-NCOa)
       Aa(row,V1+iv) = temp
       temp = 0.0D0
     enddo
     enddo
   enddo
   deallocate(B)

   ! BETA
   Nvirt = M - NCOb
   NCOc = NCOb + 1
   allocate(B(NCOb,M)); B = 0.0D0
   temp = 0.0D0
   do iv=1,Nvec
     call multlr(Cocc_transB,Fb(:,:,iv),B,NCOb,M,M,1.0D0,0.0D0)
     do i=1,NCOb
     do j=NCOc,M
       do k=1,M
         temp = temp + B(NCOc-i,k) * Cb(k,j)
       enddo
       row = (i-1) * Nvirt + (j-NCOb)
       Ab(row,V1+iv) = temp
       temp = 0.0D0
     enddo
     enddo
   enddo
   deallocate(B)
end subroutine open_MtoIANV

