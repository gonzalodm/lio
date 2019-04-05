subroutine open_RitzObtain(VMOa,VMOb,Vec,RVecA,RVecB, &
                     Ndim,Ndima,Ndimb,Sdim,Nstat)
   implicit none

   integer, intent(in) :: Ndim, Ndima, Ndimb, Sdim, Nstat
   real*8, intent(in) :: VMOa(Ndim,Sdim), VMOb(Ndim,Sdim)
   real*8, intent(in) :: Vec(Sdim,Nstat)
   real*8, intent(out) :: RVecA(Ndim,Nstat), RVecB(Ndim,Nstat)

   integer :: ii, jj, kk
   real*8 :: val
   
   RvecA = 0.0D0
   RvecB = 0.0D0

   do ii=1,Nstat
   do jj=1,Sdim
      val = Vec(jj,ii)
      ! ALPHA
      do kk=1,Ndima
         RVecA(kk,ii) = RVecA(kk,ii) + val * VMOa(kk,jj)
      enddo
      ! BETA
      do kk=1,Ndimb
         RVecB(kk,ii) = RVecB(kk,ii) + val * VMOb(kk,jj)
      enddo
   enddo
   enddo
end subroutine open_RitzObtain
