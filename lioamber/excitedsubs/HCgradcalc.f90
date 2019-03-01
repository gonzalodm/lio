subroutine HCgradcalc(rho,for,M)
use garcha_mod, only: RMM, d, r, Iz, natom, ntatom
use faint_cpu, only: int1G

   implicit none
   
   integer, intent(in) :: M
   real*8, intent(in) :: rho(M,M)
   real*8, intent(out) :: for(natom,3)

   integer :: i, j, ind, MM 
   real*8, dimension(:), allocatable :: rho_vec

   ! empaquetamos la rho total (GS + excitation)
   MM=M*(M+1)/2
   allocate(rho_vec(MM))
   ind = 1
   do i=1,M
      rho_vec(ind) = rho(i,i)
      ind = ind + 1
      do j=i+1,M
         rho_vec(ind) = rho(i,j) * 2.0D0
         ind = ind + 1
      enddo
   enddo

   call int1G(for,rho_vec,d,r,Iz,natom,ntatom)
   deallocate(rho_vec)
end subroutine HCgradcalc
