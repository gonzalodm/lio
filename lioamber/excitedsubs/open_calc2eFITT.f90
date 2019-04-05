subroutine open_calc2eFITT(Mata,Matb,Focka,Fockb,numvec,M)
use garcha_mod, only: RMM, Md, OPEN
use faint_cpu, only: int3lu
   implicit none

   integer, intent(in) :: numvec, M
   real*8, intent(inout) :: Mata(M,M,numvec), Matb(M,M,numvec)
   real*8, intent(out) :: Focka(M,M,numvec), Fockb(M,M,numvec)

   integer :: MM, MMd, M7, M9, ist, i, j
   real*8, dimension(:), allocatable :: Vec, Fock_vec

   real*8 :: notE = 0.0D0 ! This variable is not reference 
   real*8, dimension(:), allocatable :: Fmat_b, Hmat

   print*, "Calculated 2e-integrals with fitting"

   MM = M*(M+1)/2
   allocate(Vec(MM),Fock_vec(MM))

!  Hmat = this initialize Fock, in this case this is zero
   allocate(Fmat_b(MM),Hmat(MM))
   Hmat = 0.0D0
   Fmat_b = 0.0D0

!  Pointers to Gmat and Ginv
   MMd = Md * (Md+1) / 2
   M7  = 1 + 3*MM  ! now Gmat
   M9  = M7 + MMd ! now Ginv

   Mata = Mata + Matb
   do ist=1,numvec
      do i=1,M
      do j=1,i-1
         Mata(i,j,ist) = Mata(i,j,ist) + Mata(j,i,ist)
      enddo
      enddo
      call sprepack('L',M,Vec,Mata(:,:,ist))
      call int3lu(notE,Vec,Fmat_b,Fock_vec,RMM(M7:M7+MMd),&
                  RMM(M9:M9+MMd),Hmat,OPEN)
         ! int3lu(E2, rho, Fmat_b, Fmat, Gmat, 
         !        Ginv, Hmat, open_shell)

!     Transform vec -> mat
      call spunpack('L', M, Fock_vec, Focka(:,:,ist))
      call spunpack('L', M, Fmat_b, Fockb(:,:,ist))
   enddo

   deallocate(Fmat_b,Hmat,Vec,Fock_vec)
end subroutine open_calc2eFITT
