subroutine basis_deinitLR()
use lrdata, only: Coef_trans, Cocc, &
                   Cocc_trans, Cvir, Cvir_trans
   implicit none
   deallocate(Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans)
end subroutine basis_deinitLR
