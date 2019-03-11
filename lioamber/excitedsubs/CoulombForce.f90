subroutine CoulombForce(rhoExc,rhoDif,Xmat,for,M)
use garcha_mod, only: natom, Md, RMM
   implicit none

   integer, intent(in) :: M
   real*8, intent(in) :: rhoExc(M,M), rhoDif(M,M), Xmat(M,M)
   real*8, intent(inout) :: for(natom,3)

   integer :: i, j, ind, MM, M9, MMd
   real*8, dimension(:), allocatable :: af_Exc, af_Dif, af_X
   real*8, dimension(:), allocatable :: E_vec, D_vec, X_vec

   ! Pointers to Ginv
   MM = M*(M+1)/2 
   MMd = Md*(Md+1)/2
   M9 = 1+3*MM+MMd

   ! Put Densities in vector form
   allocate(E_vec(MM),D_vec(MM),X_vec(MM))
   ind = 1
   do i=1,M
      E_vec(ind) = rhoExc(i,i)
      D_vec(ind) = rhoDif(i,i) * 2.0D0
      X_vec(ind) = Xmat(i,i) * 2.0D0
      ind = ind + 1
      do j=i+1,M
         E_vec(ind) = rhoExc(i,j) * 2.0D0
         D_vec(ind) = (rhoDif(i,j) + rhoDif(j,i)) * 2.0D0 ! modifique por 2
         X_vec(ind) = (Xmat(i,j) + Xmat(j,i)) * 2.0D0 ! modifique por 2 
         ind = ind + 1
      enddo
   enddo

   ! Calculation of Variational coef with differents densities
   allocate(af_Exc(Md),af_Dif(Md),af_X(Md))
   af_Exc = 0.0D0; af_Dif = 0.0D0; af_X = 0.0D0
   call varcoef_calc(E_vec,D_vec,X_vec, & ! inputs
                     af_Exc,af_Dif,af_X, & ! outputs
                     RMM(M9:M9+MMd),MM,MMd,Md) ! Ginv and dimensions

   ! Calculate gradients of terms (P|Q)
   call int2GExc(af_Exc,af_Dif,af_X,for,natom,Md)

   ! Calculate gradients of terms (uv|P)
   call int3GExc(af_Exc,af_Dif,af_X,E_vec,D_vec,X_vec,for,natom,Md,MM)

   deallocate(af_Exc,af_Dif,af_X)
   deallocate(E_vec,D_vec,X_vec)
end subroutine CoulombForce
