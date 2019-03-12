subroutine RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,Qvec,NCO,Nvirt,Ndim)
   implicit none

   integer, intent(in) :: NCO, Nvirt, Ndim
   real*8, intent(in) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
   real*8, intent(in) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)
   real*8, intent(in) :: X(Ndim)
   real*8, intent(out) :: Rvec(Ndim), Qvec(Ndim)

   integer :: i, a, b, j, posf, pos1, NCOc
   real*8 :: temp1, temp2

   temp1 = 0.0D0; temp2 = 0.0D0
   NCOc = NCO + 1

   do i=1,NCO
   do a=1,Nvirt
      posf = (i-1) * Nvirt + a
      ! VIRTUAL PART
      do b=1,Nvirt
         pos1 = (i-1) * Nvirt + b
         temp1 = temp1 + X(pos1) * FXAB(a,b)
      enddo
      ! OCC PART
      do j=1,NCO
         pos1 = (j-1) * Nvirt + a
         temp2 = temp2 + X(pos1) * FXIJ(NCOc-i,NCOc-j)
      enddo
      Qvec(posf) = temp2
      ! OBRAIN RIA IN VECTOR FORM
      Rvec(posf) = temp2 - (temp1 + FTIA(NCOc-i,a) + 2.0D0*GXCIA(NCOc-i,a))
      temp1 = 0.0D0; temp2 = 0.0D0
   enddo
   enddo
end subroutine RCalculate
