subroutine matMOtomatAO(MatMO,MatAO,C,M,numvec,transp)
use lrdata, only: Coef_trans

   implicit none

   integer, intent(in) :: M, numvec
   logical, intent(in) :: transp
   real*8, intent(in) :: MatMO(M,M,numvec)
   real*8, intent(in) :: C(M,M)
   real*8, intent(out) :: MatAO(M,M,numvec)

   integer :: i, j, k
   real*8, dimension(:,:), allocatable :: scratch

!  CHANGE BASIS MO -> AO
   allocate(scratch(M,M))
   do i=1,numvec
      call multlr(MatMO(:,:,i),Coef_trans,scratch,M,M,M,1.0D0,0.0D0)
      call multlr(C,scratch,MatAO(:,:,i),M,M,M,1.0D0,0.0D0)
   enddo

!  WE TRANSPOSE MATRIX FOR USE IT IN C
   if ( transp ) then
     do i=1,numvec
       scratch = MatAO(:,:,i)
       do j=1,M
       do k=1,M
         MatAO(k,j,i) = scratch(j,k)
         print*, k,j,i,MatAO(k,j,i)
       enddo
       enddo
     enddo
   endif

   deallocate(scratch)
end subroutine matMOtomatAO

