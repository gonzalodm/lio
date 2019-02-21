subroutine TransDipole(Tdens,Tdip,M,Nstat)
   use garcha_mod, only: RMM
   implicit none

   integer, intent(in) :: M, Nstat
   real*8, intent(inout) :: Tdens(M,M,Nstat)
   real*8, intent(out) :: Tdip(Nstat,3)

   integer :: i, j, ist
   real*8, dimension(:,:), allocatable :: RhoRmm

!  SAVE RHO FROM RMM
   allocate(RhoRmm(M,M))
   call spunpack_rho('L',M,RMM,RhoRmm)

   do ist=1,Nstat
      do i=1,M
      do j=1,i-1
         Tdens(i,j,ist) = Tdens(i,j,ist) + Tdens(j,i,ist)
      enddo
      enddo
      call sprepack('L',M,RMM,Tdens(:,:,ist))
      call dip(Tdip(ist,:))
   enddo
   Tdip = Tdip * 2.0D0 / dsqrt(2.0D0)

!  COPY RHO OLD INTO RMM
   do i=1,M
   do j=1,i-1
      RhoRMM(i,j) = RhoRMM(i,j) * 2.0D0
   enddo
   enddo
   call sprepack('L',M,RMM,RhoRMM)
   deallocate(RhoRMM)
end subroutine TransDipole
