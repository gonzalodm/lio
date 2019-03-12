subroutine WSgradcalc(We,for,M)
use garcha_mod, only: RMM, r, d, ntatom, natom, &
                      Md
 use faint_cpu, only: intSG
   implicit none
   
   integer, intent(in) :: M
   real*8, intent(inout) :: We(M,M)
   real*8, intent(inout) :: for(natom,3)

   integer :: i, j, ind, MM, M15, MMd
   real*8, dimension(:), allocatable :: Wtot

   ! Pointers to W from excited states. This is -W
   MM = M  * (M  +1) / 2
   MMd = Md * (Md +1) / 2
   M15 = 1 + M + 4*MM + 2*MMd

   allocate(Wtot(MM))
   ind = 1
   do i=1,M
     Wtot(ind) = RMM(M15-1+ind) - We(i,i)
     ind = ind + 1
     do j=i+1,M
       Wtot(ind) = RMM(M15-1+ind) - We(i,j)*2.0D0
       ind = ind + 1
     enddo
   enddo

   ! We calculate gradients. -WS
   call intSG(for,Wtot,r,d,natom,ntatom)
   deallocate(Wtot)
end subroutine WSgradcalc
