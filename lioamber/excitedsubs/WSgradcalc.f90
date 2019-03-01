subroutine WSgradcalc(We,for,M)
use garcha_mod, only: RMM, r, d, ntatom, natom, &
                      Md
 use faint_cpu, only: intSG
   implicit none
   
   integer, intent(in) :: M
   real*8, intent(inout) :: We(M,M)
   real*8, intent(inout) :: for(natom,3)

   integer :: i, j, ind
   integer :: MM, M15, MMd
   real*8, dimension(:), allocatable :: Wtot

   print*, "calc gradiente WS total"
 
   ! pointers para sacar W del estado fundamental
   MM = M  * (M  +1) / 2
   MMd = Md * (Md +1) / 2
   M15 = 1 + M + 4*MM + 2*MMd

   print*, "W del estado fundamental, vec"
   do i=M15,M15+MM-1
      print*, RMM(i)
   enddo
   print*, "Wex mat"
   do i=1,M
   do j=1,M
      print*, i,j,We(i,j)
   enddo
   enddo

   ! Wgs ya tiene su signo menos, por eso la resta con Wexc
   allocate(Wtot(MM))
   ind = 1
!  do i=1,M
!     We(i,i) = We(i,i) * 0.5D0
!     do j=1,i
!        Wtot(ind) = RMM(M15-1+ind) - (We(i,j)*2.0D0)
!        print*, ind,RMM(M15-1+ind), We(i,j)
!        ind = ind + 1
!     enddo
!     We(i,i) = We(i,i) * 2.0D0
!  enddo
   do i=1,M
     Wtot(ind) = RMM(M15-1+ind) - We(i,i)
     print*,"d",ind,RMM(M15-1+ind), We(i,i)/2.0D0
     ind = ind + 1
     do j=i+1,M
       Wtot(ind) = RMM(M15-1+ind) - We(i,j)*2.0D0
       print*, ind,RMM(M15-1+ind)/2.0D0, We(i,j)
       ind = ind + 1
     enddo
   enddo
   !stop
   
!  print*, "Wtot = Wgs + We, vec"
!  do i=1,MM
!     print*, Wtot(i)
!  enddo

   ! calculamos el termino de grad -WS
   call intSG(for,Wtot,r,d,natom,ntatom)
end subroutine WSgradcalc
