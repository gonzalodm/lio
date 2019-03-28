module lrdata
   implicit none

   logical :: lresp = .false.
   ! Frozen Core and Valence Approximation
   logical :: FCA = .false.
   integer :: nfo = 3
   integer :: nfv = 3
   ! Change Basis
   real*8, dimension(:,:), allocatable :: Coef_trans, Cocc, Cocc_trans, &
                                          Cvir, Cvir_trans
   ! Needed for LIBINT
   real*8, dimension(:,:), allocatable :: cbas, cbasx
   logical :: fbas = .true.
   
   integer :: nstates = 3
   integer :: root = 0
   logical :: fitLR = .false.

   ! Excited States FORCES
   logical :: excited_forces = .false.
   real*8, dimension(:,:), allocatable :: forEXC


!  OPEN LINEAR RESPONSE
   real*8, dimension(:,:), allocatable :: coef_a
  
   ! Change basis
   real*8, dimension(:,:), allocatable :: Coef_transB, CoccB, Cocc_transB, &
                                          CvirB, Cvir_transB
end module lrdata