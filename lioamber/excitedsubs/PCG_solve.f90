subroutine PCG_solve(bvec,Rho_urel,Coef,E,M,NCO,Nvirt,Ndim)
use lrdata, only: cbas, fitLR
   implicit none

   integer, intent(in) :: NCO, Nvirt, Ndim, M
   real*8, intent(in) :: bvec(Ndim), E(M), Coef(M,M), Rho_urel(M,M)

   integer :: i, j, iter, maxIter
   real*8 :: beta, alpha
   logical :: conv = .false.
   real*8, dimension(:), allocatable :: R, Z, Pk, Mprec, ApIA, X
   real*8, dimension(:,:), allocatable :: Pmat, F2e, Fxc, Ftot, CopyP

!  START PRECONDITINED CONJUGATE GRADIENT
   maxIter = 50

!  INITIAL GUESS: Xo = 0
   allocate(R(Ndim))
   R = bvec

!  CALCULATE PRECONDITIONED M^(-1)
   allocate(Mprec(Ndim))
   call Prec_calculate(E,Mprec,M,NCO,Ndim)

   allocate(Pk(Ndim)); Pk = 0.0D0
   call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

   allocate(Pmat(M,M),F2e(M,M),Fxc(M,M),Ftot(M,M))
   allocate(ApIA(Ndim),X(Ndim),CopyP(M,M)); X = 0.0D0

   write(*,*) ""
   write(*,"(1X,A)") "Start PCG loop"
   write(*,*) ""
   do iter=1,maxIter

!  CONVERT TRIAL VECTORS TO AO BASIS
   call VecToMat(Pk,Pmat,Coef,Ndim,NCO,M)

!  CALCULATE TWO ELECTRON PART
   F2e = 0.0D0
   if ( .not. fitLR ) then
      call g2g_calculate2e(Pmat,cbas,1,F2e,1)
      F2e = 2.0D0 * F2e
   else
      CopyP = Pmat
      call calc2eFITT(CopyP,F2e,1,M)
   endif

!  CALCULATE XC PART
   Fxc = 0.0D0
   call g2g_calculateg(Pmat,Fxc,2)

!  OBTAIN FOCK TOTAL AND ADD TERM (Ea-Ei)Pk AND
!  CHANGE BASIS AO -> MO
   call total_fock(F2e,Fxc,Ftot,M)
   call Ap_calculate(Ftot,Pk,Coef,E,ApIA,M,NCO,Nvirt,Ndim)

!  CALCULATE ALPHA
   call Alpha_calc(Pk,ApIA,alpha,Ndim)

!  NEW X VECTOR
   X = X + alpha * Pk

!  NEW R VECTOR
   R = R - alpha * APIA

!  CHECK CONVERGENCE
   call error(R,conv,Ndim,iter)
   if (conv .eqv. .true.) then
      write(*,"(1X,A,1X,I2,1X,A)") "Convergence achieved in",iter,"itertaions"
      exit
   endif

!  GET NEW BETA AND Pk
   call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

   enddo ! ENDDO LOOP PCG

!  FORM RELAXED DENSITY OF EXCITED STATE
   call RelaxedDensity(X,Rho_urel,Coef,M,NCO,Ndim)

   deallocate(R,Mprec,Pk,Pmat,F2e,Fxc,Ftot,ApIA,X,CopyP)
end subroutine PCG_solve
