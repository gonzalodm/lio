subroutine open_subspaceMat(Aa,Ab,Va,Vb,H,NdimA,NdimB,Ndim,Sdim)
   implicit none

   integer, intent(in) :: NdimA, NdimB, Ndim, Sdim
   real*8, intent(in) :: Aa(Ndim,Sdim), Ab(Ndim,Sdim)
   real*8, intent(in) :: Va(Ndim,Sdim), Vb(Ndim,Sdim)
   real*8, intent(out) :: H(Sdim,Sdim)

   integer :: ii, jj
   real*8 :: temp_a, temp_b

   do ii=1,Sdim
   do jj=1,Sdim
      call prod(Va(:,ii),Aa(:,jj),temp_a,Ndima)
      call prod(Vb(:,ii),Ab(:,jj),temp_b,Ndimb)
      H(ii,jj) = temp_a + temp_b
   enddo
   enddo
end subroutine open_subspaceMat

subroutine prod(A,B,temp,N)
   implicit none
 
   integer, intent(in) :: N
   real*8, intent(in) :: A(N), B(N)
   real*8, intent(out) :: temp

   integer :: ii
   
   temp = 0.0D0
   
   do ii=1,N
      temp = temp + A(ii) * B(ii)
   enddo
end subroutine prod
