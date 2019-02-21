module lrdata
!Nvirt = number of virtual molecular orbitals
!dim = dimension of matrix Linear Response (Nvirt x NCO)
!cbas,cbasx = Needed for libint (Temporary)
!nstates = number of excited states to calculate
!root = excited state chosen for optimization
   implicit none

   logical :: lresp = .false.
   logical :: fbas = .true.
   logical :: fitLR = .false.
   integer :: Nvirt, dim, NCOlr, Mlr
   integer :: nstates = 3
   integer :: root = 0
   ! PASOS EN LR EN DINAMICA
   integer :: StepLR = 0
   real*8, dimension(:,:), allocatable :: eigenvec, cbas, cbasx
   real*8, dimension(:), allocatable :: eigenval
! Use Frozen Core Approximation
! nfo = number of molecular orbitals occupied with lowest energy deleted
! nfv = number of molecular orbiatls virtual with higher energy deleted
   logical :: FCA = .false.
   integer :: nfo = 3
   integer :: nfv = 3
! Matrix needed for change basis
   real*8, dimension(:,:), allocatable :: Coef_trans, Cocc, Cocc_trans, &
                                          Cvir, Cvir_trans
end module lrdata
