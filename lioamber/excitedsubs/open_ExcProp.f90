subroutine open_ExcProp(Calpha,Cbeta,Ealpha,Ebeta,NCOa,NCOb)
use garcha_mod, only: M
   implicit none
   
   integer, intent(inout) :: NCOa, NCOb
   real*8, intent(inout) :: Calpha(:,:), Cbeta(:,:)
   real*8, intent(inout) :: Ealpha(:), Ebeta(:)
   
   integer :: i, j, nada


   ! DEBUG GAMES, en games NCOa = 2, NCOb = 1
   ! EN LIO es al reves NCOa =1, NCOb =2
   open(unit=456,file="coefab")
   do i=1,M
   do j=1,M
      read(456,*) nada, Cbeta(j,i), Calpha(j,i)
   enddo
   enddo
   close(456)
   open(unit=456,file="eneab")
   do i=1,M
      read(456,*) nada, Ebeta(i), Ealpha(i)
   enddo
   close(456)

   ! BASIS INITIALIZATION
   call open_basis_initLR(Calpha,Cbeta,NCOa,NCOb,M)

   call open_linear_response(Calpha,Cbeta,Ealpha,Ebeta, &
                              M,NCOa,NCOb)
  
!  print*, "coefientes alpha y beta"
!  do i=1,M
!     print*, "ener a b", i, Ealpha(i),Ebeta(i)
!  do j=1,M
!     print*, j,Calpha(i,j),Cbeta(i,j)
!  enddo
!  enddo

   call open_basis_deinitLR()
   


end subroutine open_ExcProp
