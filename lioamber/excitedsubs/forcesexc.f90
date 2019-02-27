subroutine forcesexc(Dgs,DiffExc,Xmat,Cof,M,NCO)
use garcha_mod, only: natom
! Dgs = matriz densidad del estado fundamental
! DiffExc = matriz densidad diferencia relajada
! Xmat = vectores Linear Response en AO normalizados a 1/2
! Cof = coeficientes de los orbitales moleculares
   implicit none
  
   integer, intent(in) :: M, NCO
   real*8, intent(in) :: Dgs(M,M), DiffExc(M,M), Xmat(M,M), Cof(M,M)

   integer :: i, j
   real*8, dimension(:,:), allocatable :: forces_exc
   print*, "Calculate excited state gradients"
   do i=1,M
   do j=1,M
      !print*, i,j,Dgs(i,j),DiffExc(i,j),Cof(i,j)
      print*, i,j,Cof(i,j)
   enddo
   enddo

!  print*, ""
!  do i=1,M
!  do j=1,M
!     print*, i,j,Xmat(i,j)
!  enddo
!  enddo
   
   allocate(forces_exc(3,natom))
   call g2g_calcgrdexc(DiffExc,Xmat,Cof,forces_exc,NCO)
   deallocate(forces_exc)

   


end subroutine forcesexc
