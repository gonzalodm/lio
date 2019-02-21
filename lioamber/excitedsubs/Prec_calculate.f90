subroutine Prec_calculate(Ener,Mprec,M,NCO,Ndim)
   implicit none

   integer, intent(in) :: M, NCO, Ndim
   real*8, intent(in) :: Ener(M)
   real*8, intent(out) :: Mprec(Ndim)

   integer :: i, j, Nvirt, NCOc, pos
   real*8 :: temp

   Nvirt = M - NCO
   NCOc = NCO + 1

   do i=1,NCO
   do j=1,Nvirt
      pos = (i-1) * Nvirt + j
      temp = Ener(j+NCO) - Ener(NCOc-i)
      temp = 1.0D0 / temp
      Mprec(pos) = temp
   enddo
   enddo
end subroutine Prec_calculate
