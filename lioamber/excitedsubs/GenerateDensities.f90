subroutine GenerateDensities(X,C,P,T,M,Ndim,NCO,Nvirt)
   implicit none

   integer, intent(in) :: M, Ndim, NCO, Nvirt
   real*8, intent(in) :: X(Ndim), C(M,M)
   real*8, intent(out) :: P(M,M), T(M,M)

   real*8, dimension(:,:), allocatable :: PMO

   ! CALCULATE DIFFERENCE UNRELAXED DENSITY MATRIX
   allocate(PMO(M,M))
   call UnDiffDens(X,PMO,NCO,Nvirt,M,Ndim)
   call matMOtoMatAO(PMO,P,C,M,1,.false.)
   deallocate(PMO)

   ! CALCULATE TRANSITION DENSITY MATRIX
   call XmatForm(X,C,T,Ndim,NCO,Nvirt,M)
end subroutine GenerateDensities
