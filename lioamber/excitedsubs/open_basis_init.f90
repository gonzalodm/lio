subroutine open_basis_initLR(Ca,Cb,NCOa,NCOb,M)
!  This subroutine initializes the matrix needed for 
!  the change basis in LR, Zvector, Forces calculations.
use lrdata, only: Coef_trans, Cocc, Cocc_trans, & ! alpha
                  Cvir, Cvir_trans, &
                  Coef_transB, CoccB, Cocc_transB, & ! beta
                  CvirB, Cvir_transB
   implicit none
   integer, intent(in) :: M, NCOa, NCOb
   real*8, intent(in) :: Ca(M,M), Cb(M,M)

   integer :: Nvirta, Nvirtb, j

   Nvirta = M - NCOa
   Nvirtb = M - NCOb

   ! allocate dimension alpha
   allocate(Coef_trans(M,M),Cocc(M,NCOa),Cocc_trans(NCOa,M), &
            Cvir(M,Nvirta),Cvir_trans(Nvirta,M))

   ! allocate dimension beta
   allocate(Coef_transB(M,M),CoccB(M,NCOb),Cocc_transB(NCOb,M), &
            CvirB(M,Nvirtb),Cvir_transB(Nvirtb,M))

   ! ALPHA
   do j=1,NCOa
      Cocc(:,j) = Ca(:,j)
      Cocc_trans(j,:) = Ca(:,j)
   enddo
   do j=1,Nvirta
      Cvir(:,j) = Ca(:,NCOa+j)
      Cvir_trans(j,:) = Ca(:,NCOa+j)
   enddo
   Coef_trans = transpose(Ca)

   ! BETA
   do j=1,NCOb
      CoccB(:,j) = Cb(:,j)
      Cocc_transB(j,:) = Cb(:,j)
   enddo
   do j=1,Nvirtb
      CvirB(:,j) = Cb(:,NCOb+j)
      Cvir_transB(j,:) = Cb(:,NCOb+j)
   enddo
   Coef_transB = transpose(Cb)
end subroutine open_basis_initLR
