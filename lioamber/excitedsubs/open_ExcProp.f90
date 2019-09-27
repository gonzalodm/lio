subroutine open_ExcProp(Calpha,Cbeta,Ealpha,Ebeta,NCOa,NCOb)
use lrdata, only: second_LR, doing_SLR
use garcha_mod, only: M
   implicit none
   
   integer, intent(inout) :: NCOa, NCOb
   real*8, intent(inout) :: Calpha(:,:), Cbeta(:,:)
   real*8, intent(inout) :: Ealpha(:), Ebeta(:)
   
   integer :: i, j, nada

   ! BASIS INITIALIZATION
   call open_basis_initLR(Calpha,Cbeta,NCOa,NCOb,M)

   if (second_LR) doing_SLR = .false.

   ! TDA CALCULATION
   call open_linear_response(Calpha,Cbeta,Ealpha,Ebeta, &
                              M,NCOa,NCOb)

   if (second_LR) then
      second_LR = .false.
      doing_SLR = .true.
      call open_second_LinearResponse(Calpha,Cbeta,Ealpha,Ebeta,&
                                      M,NCOa,NCOb)
   endif

   ! BASIS DEINITIALIZATION
   call open_basis_deinitLR()

end subroutine open_ExcProp
