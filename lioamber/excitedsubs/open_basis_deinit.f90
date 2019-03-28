subroutine open_basis_deinitLR()
use lrdata, only: Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans, &
                  Coef_transB, CoccB, Cocc_transB, CvirB, Cvir_transB

   implicit none
   deallocate(Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans)
   deallocate(Coef_transB, CoccB, Cocc_transB, CvirB, Cvir_transB)
end subroutine open_basis_deinitLR
