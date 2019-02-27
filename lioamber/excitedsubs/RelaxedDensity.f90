subroutine RelaxedDensity(Z,Rho_urel,Xmat,C,M,NCO,N)
use garcha_mod, only: RMM, &
                      Smat, mulliken, Iz, natom, Nuc ! for mulliken
use fileio    , only: write_population ! for mulliken
use lrdata    , only: root ! for mulliken

   implicit none

   integer, intent(in) :: M, NCO, N
   real*8, intent(in) :: Z(N), C(M,M), Rho_urel(M,M), Xmat(M,M)

   integer :: i, j, Nvirt, NCOc, pos, M2
   real*8, dimension(:,:), allocatable :: Rho_fund, Rel_diff, Rho_exc
   real*8, dimension(:,:), allocatable :: Zmo, Zao, RhoRMM
   real*8, dimension(:), allocatable :: qv ! for mulliken

!  EXTRACT RHO FUND FROM RMM
   allocate(Rho_fund(M,M),RhoRMM(M,M))
   call spunpack_rho('L',M,RMM,Rho_fund)

!  CONVERT Z IN AO BASIS     
   allocate(Zmo(M,M)); Zmo = 0.0D0
   Nvirt = M - NCO
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
   allocate(Rel_diff(M,M))
   Rel_diff = Rho_urel + Zao

   allocate(Rho_exc(M,M))
   do i=1,M
   do j=1,M
      Rho_exc(i,j) = Rho_fund(i,j) + Rel_diff(i,j) + Rel_diff(j,i)
   enddo
   enddo

   print*, "Rho fundamental"
   do i=1,M
   do j=1,M
      print*, i,j,Rho_fund(i,j)
   enddo
   enddo

   print*, "relaxed density excited state"
   do i=1,M
   do j=1,M
      print*, i,j,Rho_exc(i,j)
   enddo
   enddo

!  SAVE EXCITED RHO INTO RMM
   if (.false.) then
      M2 = M * 2
      do i=1,M
         RMM(i + (M2-i)*(i-1)/2) = Rho_exc(i,i)
         do j = i+1, M
            RMM(j + (M2-i)*(i-1)/2) = Rho_exc(i,j) * 2.0D0
         enddo
      enddo
   endif

! Temporary for forces in excited states
   call forcesexc(Rho_fund,Rel_diff,Xmat,C,M,NCO)
! ######################################

!   TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
!   MULLIKEN POPULATION
!   Hay un bug con mulliken, esto imprime las cargas de mulliken
!   del estado excitado root, pero debido a que guardo en RMM la
!   matriz de estado excitado el archivo mulliken tambien tiene
!   las cargas de mull del estado excitado
!   TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
   if (mulliken) then
      open(unit=456,file="mulliken.exc")
      write(456,"(A,1X,I2)") "Mulliken population from excited state", root
      allocate(qv(natom))
      do i=1,natom
         qv(i) = real(Iz(i))
      enddo
      call mulliken_calc(natom, M, Rho_exc, Smat, Nuc, qv)
      call write_population(natom, Iz, qv, 0, 456)
      deallocate(qv)
      close(456)
   endif

   deallocate(Rho_exc,Rel_diff)
end subroutine RelaxedDensity
