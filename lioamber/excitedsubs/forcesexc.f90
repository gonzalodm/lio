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
   real*8, dimension(:,:), allocatable :: forXC, forWS, forHV, forCou, forTot
   real*8, dimension(:,:), allocatable :: Wexc, rhoTot

   print*, "Calculate excited state gradients"

   Ndim = NCO * (M - NCO)
   ! Calculate gradients of overlap
   allocate(Wexc(M,M)); Wexc = 0.0D0
   call Wcalculate(Zvec,DiffExc,Cof,Wexc,Ndim,M)
   allocate(forWS(natom,3)); forWS = 0.0D0
   call WSgradcalc(Wexc,forWS,M)
   deallocate(Wexc)

   ! Calculate total rho of excited state
   allocate(rhoTot(M,M))
   rhoTot = Dgs + DiffExc + transpose(DiffExc)

   ! Calculate gradients of Core and Nuclear
   allocate(forHV(natom,3)); forHV = 0.0D0
   call HVgradcalc(rhoTot,forHV,M)

   ! Calculate gradients of Coulomb
   allocate(forCou(natom,3)); forCou = 0.0D0
   call CoulombForce(rhoTot,DiffExc,Xmat,forCou,M)
   deallocate(rhoTot)
   
   ! Forces XC in excited state total ( DFT + excitation )
   allocate(forXC(3,natom)); forXC = 0.0D0
   call g2g_calcgrdexc(DiffExc,Xmat,Cof,forXC,NCO)

   allocate(forTot(natom,3))
   forTot = forWS + forHV + forCou + transpose(forXC)

   print*, "TOTAL FORCES IN EXCITED STATE"
   do i=1,natom
      print*, i, forTot(i,1), forTot(i,2), forTot(i,3)
   enddo

   deallocate(forXC,forWS,forHV,forCou,forTot)
end subroutine forcesexc
