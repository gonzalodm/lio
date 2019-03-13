subroutine Zvector(C,Ene,X,TundAO,Xmat,Zvec,Qvec,Gxc,NCO,M,Ndim,Nvirt)
use lrdata, only: cbas, root, fitLR
   implicit none

   integer, intent(in) :: NCO, M, Ndim, Nvirt
   real*8, intent(in) :: C(M,M), Ene(M), TundAO(M,M)
   real*8, intent(in) :: Xmat(M,M), X(Ndim)
   real*8, intent(inout) :: Zvec(Ndim), Qvec(Ndim), Gxc(M,M)

   real*8, dimension(:,:), allocatable :: FX, FT
   real*8, dimension(:,:,:), allocatable :: F2e, PA
   real*8, dimension(:,:), allocatable :: FXAB, FXIJ, FTIA, GXCIA
   real*8, dimension(:), allocatable :: Rvec

   print*, ""
   print*,"======================================="
   print*,"            Z-VECTOR METHOD"
   print*,"======================================="
   print*, ""

   write(*,"(1X,A,1X,I2)") "FORM RELAXED DENSITY MATRIX FOR EXCITED STATE:", root
   write(*,*) ""

!  CALCULATE THIRD DERIVATIVE FOCK
   Gxc = 0.0D0
   call g2g_calculateg(Xmat,Gxc,3)

!  CALCULATE SECOND DERIVATIVE FOCK
   allocate(FX(M,M)); FX = 0.0D0
   call g2g_calculateg(Xmat,FX,2)
   allocate(FT(M,M)); FT = 0.0D0
   call g2g_calculateg(TundAO,FT,2)

!  CALCULATE TWO ELECTRON FOCK
!  1 = tiene la parte 2e de Xmat
!  2 = tiene la parte 2e de Tund
   allocate(PA(M,M,2),F2e(M,M,2)); F2e = 0.0D0
   PA(:,:,1) = Xmat; PA(:,:,2) = TundAO

   if ( .not. fitLR) then
      call g2g_calculate2e(PA,cbas,2,F2e,1)
      F2e = 2.0D0 * F2e
   else
      call calc2eFITT(PA,F2e,2,M)
   endif
   deallocate(PA)

!  FT + F2eT
   FT = F2e(:,:,2) + 2.0D0 * FT
!  FX + F2eX
   FX = F2e(:,:,1) + 2.0D0 * FX
   deallocate(F2e)

!  CHANGE BASIS OF ALL FOCK TYPE MATRIX
   allocate(FXAB(Nvirt,Nvirt),FXIJ(NCO,NCO))
   allocate(FTIA(NCO,Nvirt),GXCIA(NCO,Nvirt))
   call ChangeBasisF(FX,FT,Gxc,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
   deallocate(FX,FT)

!  CALCULATE VECTOR R (A * X = R)
   allocate(Rvec(Ndim))
   call RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,Qvec,NCO,Nvirt,Ndim)

!  SOLVE EQUATION AX=R WITH PCG METHOD     
   call PCG_solve(Rvec,C,Ene,Zvec,M,NCO,Nvirt,Ndim)

   deallocate(FXAB,FXIJ,FTIA,GXCIA,Rvec)
end subroutine Zvector
