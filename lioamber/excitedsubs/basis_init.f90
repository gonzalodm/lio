subroutine basis_initLR(Coef,M,NCO,Nvirt)
!  This subroutine initializes the matrix needed for 
!  the change of basis in linear response and CPKS calculations.
use lrdata, only: Coef_trans, Cocc, &
                        Cocc_trans, Cvir, Cvir_trans
   implicit none
   integer, intent(in) :: M, NCO, Nvirt
   real*8, intent(in) :: Coef(M,M)

   integer :: i, j

   allocate(Coef_trans(M,M),Cocc(M,NCO),Cocc_trans(NCO,M), &
            Cvir(M,Nvirt),Cvir_trans(Nvirt,M))

   do j=1,NCO
      Cocc(:,j) = Coef(:,j)
      Cocc_trans(j,:) = Coef(:,j)
   enddo

   do j=1,Nvirt
      Cvir(:,j) = Coef(:,NCO+j)
      Cvir_trans(j,:) = Coef(:,NCO+j)
   enddo

   Coef_trans = transpose(Coef)
end subroutine basis_initLR
