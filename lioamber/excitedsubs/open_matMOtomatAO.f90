subroutine open_matMOtomatAO(MOa,MOb,AOa,AOb,Ca,Cb,M,numvec,transp)
use lrdata, only: Coef_trans, Coef_transB

   implicit none

   integer, intent(in) :: M, numvec
   logical, intent(in) :: transp
   real*8, intent(in) :: Ca(M,M), Cb(M,M)
   real*8, intent(in) :: MOa(M,M,numvec), MOb(M,M,numvec)
   real*8, intent(out) :: AOa(M,M,numvec), AOb(M,M,numvec)

   integer :: i, j, k
   real*8, dimension(:,:), allocatable :: scratchA, scratchB

!  CHANGE BASIS MO -> AO
   allocate(scratchA(M,M),scratchB(M,M))
   do i=1,numvec
      ! alpha
      call multlr(MOa(:,:,i),Coef_trans,scratchA,M,M,M,1.0D0,0.0D0)
      call multlr(Ca,scratchA,AOa(:,:,i),M,M,M,1.0D0,0.0D0)
      ! beta
      call multlr(MOb(:,:,i),Coef_transB,scratchB,M,M,M,1.0D0,0.0D0)
      call multlr(Cb,scratchB,AOb(:,:,i),M,M,M,1.0D0,0.0D0)
   enddo

!  WE TRANSPOSE MATRIX FOR USE IT IN C
   if ( transp ) then
     do i=1,numvec
       scratchA = AOa(:,:,i)
       scratchB = AOb(:,:,i)
       do j=1,M
       do k=1,M
         AOa(k,j,i) = scratchA(j,k)
         AOb(k,j,i) = scratchB(j,k)
       enddo
       enddo
     enddo
   endif

   deallocate(scratchA,scratchB)
end subroutine open_matMOtomatAO
