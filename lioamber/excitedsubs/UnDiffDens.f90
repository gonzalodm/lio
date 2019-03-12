subroutine UnDiffDens(X,T,NCO,Nvirt,M,Ndim)
   implicit none

   integer, intent(in) :: NCO, Nvirt, M, Ndim
   real*8, intent(in) :: X(Ndim)
   real*8, intent(out) :: T(M,M)

   integer :: i, j, pos, NCOc
   real*8, dimension(:,:), allocatable :: XM, XMtrans, Ptrash

   allocate(XM(NCO,Nvirt),XMtrans(Nvirt,NCO))

!  WE NORMALIZE TO 1/2
   T = 0.0D0
   NCOc = NCO + 1

   do i=1,NCO
   do j=1,Nvirt
     pos = (i-1) * Nvirt + j
     XM(i,j) = X(pos)
     XMtrans(j,i) = X(pos)
   enddo
   enddo

!  FORM UNRELAXED DIFFERENCE DENSITY MATRIX
   allocate(Ptrash(NCO,NCO))
!  FORM BLOCK OCC-OCC
   call multlr(XM,XMtrans,Ptrash,NCO,Nvirt,NCO,-1.0D0,0.0D0)
   do i=1,NCO
   do j=i,NCO
      T(NCOc-i,NCOc-j) = Ptrash(i,j)
      T(NCOc-j,NCOc-i) = Ptrash(i,j)
   enddo
   enddo

!  FORM BLOCK VIR - VIR
   deallocate(Ptrash); allocate(Ptrash(Nvirt,Nvirt))
   call multlr(XMtrans,XM,Ptrash,Nvirt,NCO,Nvirt,1.0D0,0.0D0)
   do i=1,Nvirt
   do j=i,Nvirt
      T(NCO+i,NCO+j) = Ptrash(i,j)
      T(NCO+j,NCO+i) = Ptrash(i,j)
   enddo
   enddo
   deallocate(Ptrash,XM,XMtrans)
end subroutine UnDiffDens
