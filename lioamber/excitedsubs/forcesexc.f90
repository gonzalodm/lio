subroutine forcesexc(rhoTot,DiffExc,Zvec,Xmat,&
                     Qvec,GxcAO,Xlr,Cscf,Escf,dE,&
                     forTot,M,Ndim,natom,NCO)
   implicit none
  
   integer, intent(in) :: M, Ndim, natom, NCO
   real*8, intent(in) :: rhoTot(M,M), DiffExc(M,M), Zvec(Ndim), Xmat(M,M)
   real*8, intent(in) :: Qvec(Ndim), GxcAO(M,M), Xlr(Ndim), Cscf(M,M), Escf(M)
   real*8, intent(in) :: dE
   real*8, intent(out) :: forTot(natom,3)

   integer :: i
   real*8, dimension(:,:), allocatable :: Wexc
   real*8, dimension(:,:), allocatable :: forXC, forWS, forHV, forCou

   print*, "Calculate excited state gradients"

   ! Calculate gradients of overlap
   allocate(Wexc(M,M)); Wexc = 0.0D0
   call Wcalculate(Zvec,DiffExc,Qvec,GxcAO,Xlr,Cscf,dE,Escf,Wexc,Ndim,M,NCO)
   allocate(forWS(natom,3)); forWS = 0.0D0
   call WSgradcalc(Wexc,forWS,M)
   deallocate(Wexc)

   ! Calculate gradients of Core and Nuclear
   allocate(forHV(natom,3)); forHV = 0.0D0
   call HVgradcalc(rhoTot,forHV,M)

   ! Calculate gradients of Coulomb
   allocate(forCou(natom,3)); forCou = 0.0D0
   call CoulombForce(rhoTot,DiffExc,Xmat,forCou,M)
   
   ! Forces XC in excited state total ( DFT + excitation )
   allocate(forXC(3,natom)); forXC = 0.0D0
   call g2g_calcgrdexc(DiffExc,Xmat,Cscf,forXC,NCO)

   forTot = forWS + forHV + forCou + transpose(forXC)

!  print*, "TOTAL FORCES IN EXCITED STATE"
!  do i=1,natom
!     print*, i, forTot(i,1), forTot(i,2), forTot(i,3)
!  enddo

   deallocate(forXC,forWS,forHV,forCou)
end subroutine forcesexc
