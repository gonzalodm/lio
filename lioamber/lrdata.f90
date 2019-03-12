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

   real*8, dimension(:,:), allocatable :: forEXC
end module lrdata
