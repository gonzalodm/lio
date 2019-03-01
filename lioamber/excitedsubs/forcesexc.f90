subroutine forcesexc(Dgs,DiffExc,Zvec,Xmat,Cof,M,NCO)
use garcha_mod, only: natom
! Dgs = matriz densidad del estado fundamental
! DiffExc = matriz densidad diferencia relajada
! Xmat = vectores Linear Response en AO normalizados a 1/2
! Cof = coeficientes de los orbitales moleculares
   implicit none
  
   integer, intent(in) :: M, NCO
   real*8, intent(in) :: Dgs(M,M), DiffExc(M,M), Xmat(M,M), Cof(M,M)
   real*8, intent(in) :: Zvec(NCO*(M-NCO))

   integer :: i, j, Ndim
   real*8, dimension(:,:), allocatable :: forces_exc, forWS, forHC, forCou, forTot
   real*8, dimension(:,:), allocatable :: Wexc, rhoTot
   print*, "Calculate excited state gradients"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Dgs(i,j),DiffExc(i,j),Xmat(i,j)
!     !print*, i,j,Cof(i,j)
!  enddo
!  enddo

!  print*, ""
!  do i=1,M
!  do j=1,M
!     print*, i,j,Xmat(i,j)
!  enddo
!  enddo
   Ndim = NCO * (M - NCO)
   allocate(Wexc(M,M)); Wexc = 0.0D0
   ! Wexc sale sin el signo menos
   call Wcalculate(Zvec,DiffExc,Cof,Wexc,Ndim,M)
!  print*, "W in AO"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Wexc(i,j)
!  enddo
!  enddo
   allocate(forWS(natom,3)); forWS = 0.0D0
   call WSgradcalc(Wexc,forWS,M)
   print*, "forces WS"
   do i=1,natom
      print*, i, forWS(i,1),forWS(i,2),forWS(i,3)
   enddo

   ! calculamos los gradientes del core y el nucleo
   allocate(forHC(natom,3)); forHC = 0.0D0
   allocate(rhoTot(M,M)); rhoTot = Dgs + DiffExc + transpose(DiffExc)
   call HCgradcalc(rhoTot,forHC,M)
   print*, "forces H+V"
   do i=1,natom
      print*, i,forHC(i,1),forHC(i,2),forHC(i,3)
   enddo

   allocate(forCou(natom,3)); forCou = 0.0D0
   call CoulombForce(rhoTot,DiffExc,Xmat,forCou,M)
   
   ! Forces XC in excited state total ( DFT + excitation )
   allocate(forces_exc(3,natom)); forces_exc = 0.0D0
   call g2g_calcgrdexc(DiffExc,Xmat,Cof,forces_exc,NCO)
   ! transponer forces_exc ya que las fuerzas en realidad son (natom,3)
   print*, "FORCES EXC"
   do i=1,natom
      print*, i, forces_exc(1,i), forces_exc(2,i), forces_exc(3,i)
   enddo
   allocate(forTot(natom,3))
   forTot = forWS + forHC + forCou + transpose(forces_exc)
   print*, "FORCES TOTAL"
   do i=1,natom
      print*, i, forTot(i,1), forTot(i,2), forTot(i,3)
   enddo


   deallocate(forces_exc)

   


end subroutine forcesexc
