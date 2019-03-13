subroutine RelaxedDensity(Z,Rho_urel,C,Rho_exc,Rel_diff,M,NCO,N,Nvirt)
use garcha_mod, only: RMM

   implicit none

   integer, intent(in) :: M, NCO, N, Nvirt
   real*8, intent(in) :: Z(N), C(M,M), Rho_urel(M,M)
   real*8, intent(out) :: Rho_exc(M,M), Rel_diff(M,M)

   integer :: i, j, NCOc, pos
   real*8, dimension(:,:), allocatable :: Zmo, Zao, Rho_fund

!  EXTRACT RHO FUND FROM RMM
   allocate(Rho_fund(M,M)); Rho_fund = 0.0D0
   call spunpack_rho('L',M,RMM,Rho_fund)

!  CONVERT Z IN AO BASIS     
   allocate(Zmo(M,M)); Zmo = 0.0D0
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,Nvirt
      pos = (i - 1) * Nvirt + j
      Zmo(NCOc-i,NCO+j) = Z(pos)
   enddo
   enddo
   allocate(Zao(M,M))
   call matMOtomatAO(Zmo,Zao,C,M,1,.false.)
   deallocate(Zmo)

   Rel_diff = Rho_urel + Zao
   deallocate(Zao)

   Rho_exc = Rho_fund + Rel_diff + transpose(Rel_diff)

!  print*, "Rho fundamental"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Rho_fund(i,j)
!  enddo
!  enddo
   deallocate(Rho_fund)

!  print*, "relaxed density excited state"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Rho_exc(i,j)
!  enddo
!  enddo
end subroutine RelaxedDensity
