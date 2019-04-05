subroutine open_ExcProp(Calpha,Cbeta,Ealpha,Ebeta,NCOa,NCOb)
use garcha_mod, only: M
   implicit none
   
   integer, intent(inout) :: NCOa, NCOb
   real*8, intent(inout) :: Calpha(:,:), Cbeta(:,:)
   real*8, intent(inout) :: Ealpha(:), Ebeta(:)
   
   integer :: i, j, nada

   ! BASIS INITIALIZATION
   call open_basis_initLR(Calpha,Cbeta,NCOa,NCOb,M)

   ! TDA CALCULATION
   call open_linear_response(Calpha,Cbeta,Ealpha,Ebeta, &
                              M,NCOa,NCOb)

   ! BASIS DEINITIALIZATION
   call open_basis_deinitLR()

end subroutine open_ExcProp
