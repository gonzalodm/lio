subroutine varcoef_calc(Evec,Dvec,Xvec,afE,afD,afX,Ginv,MM,MMd,Md)
use garcha_mod, only: cool, cools, kkind, kkinds, kknumd, kknums
! cool: precalculated 2e terms in double precision. (uv|P)                 !
! kkind: precalculated indexes for double precision Fock matrix elements.  !
! kknumd: number of precalculated double precision Fock matrix elements.   !

! cools: precalculated 2e terms in single precision. (uv|P)                !
! kkinds: precalculated indexes for single precision Fock matrix elements. !
! kknums: number of precalculated single precision Fock matrix elements.   !
   implicit none

   integer, intent(in) :: MM, MMd, Md
   real*8, intent(in) :: Evec(MM), Dvec(MM), Xvec(MM), Ginv(MMd)
   real*8, intent(out) :: afE(Md), afD(Md), afX(Md)

   integer :: ind, kk_ind, iikk, m_ind, k_ind
   real :: int3s
   real*8 :: int3d
   real*8, dimension(:), allocatable :: RcE, RcD, RcX

   allocate(RcE(Md),RcD(Md),RcX(Md))
   RcE = 0.0D0; RcD = 0.0D0; RcX = 0.0D0

   ! accumulate Rc double precision
   do kk_ind = 1, kknumd
      iikk = (kk_ind - 1) * Md
      do k_ind = 1, Md
         int3d = cool(iikk + k_ind)
         ind = kkind(kk_ind)
         RcE(k_ind) = RcE(k_ind) + Evec(ind) * int3d
         RcD(k_ind) = RcD(k_ind) + Dvec(ind) * int3d
         RcX(k_ind) = RcX(k_ind) + Xvec(ind) * int3d
      enddo
   enddo

   ! accumulate Rc simple precision
   do kk_ind = 1, kknums
      iikk = (kk_ind - 1) * Md
      do k_ind = 1, Md
         int3s = cools(iikk + k_ind)
         ind = kkinds(kk_ind)
         RcE(k_ind) = RcE(k_ind) + Evec(ind) * int3s
         RcD(k_ind) = RcD(k_ind) + Dvec(ind) * int3s
         RcX(k_ind) = RcX(k_ind) + Xvec(ind) * int3s
      enddo
   enddo

   ! Calculation of variational coefficients with differents densities
   do m_ind = 1, Md
      afE(m_ind) = 0.0D0
      afD(m_ind) = 0.0D0
      afX(m_ind) = 0.0D0
      do k_ind = 1, m_ind-1
         ! total Excited Density
         afE(m_ind) = afE(m_ind) + &
                  RcE(k_ind) * Ginv(m_ind + (2*Md-k_ind)*(k_ind-1)/2)

         ! Difference relaxed Density
         afD(m_ind) = afD(m_ind) + &
                  RcD(k_ind) * Ginv(m_ind + (2*Md-k_ind)*(k_ind-1)/2)

         ! Transition Density ( LR vectors )
         afX(m_ind) = afX(m_ind) + &
                  RcX(k_ind) * Ginv(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
      enddo
      do k_ind = m_ind, Md
         ! total Excited Density
         afE(m_ind) = afE(m_ind) + &
                  RcE(k_ind) * Ginv(k_ind + (2*Md-m_ind)*(m_ind-1)/2)

         ! Difference relaxed Density
         afD(m_ind) = afD(m_ind) + &
                  RcD(k_ind) * Ginv(k_ind + (2*Md-m_ind)*(m_ind-1)/2)

         ! Transition Density ( vectors of LR )
         afX(m_ind) = afX(m_ind) + &
                  RcX(k_ind) * Ginv(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
      enddo
   enddo

   deallocate(RcE,RcD,RcX)
end subroutine varcoef_calc
