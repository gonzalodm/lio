subroutine open_residual(Val,Vec,Ra,Rb,Aa,Ab,Wa,Wb,&
                         Ndima,Ndimb,Ndim,Sdim,Nstat)
   implicit none

   integer, intent(in) :: Ndima, Ndimb, Ndim, Sdim, Nstat
   real*8, intent(in) :: Val(Nstat),Vec(Sdim,Nstat)
   real*8, intent(in) :: Ra(Ndim,Nstat), Rb(Ndim,Nstat)
   real*8, intent(in) :: Aa(Ndim,Sdim), Ab(Ndim,Sdim)
   real*8, intent(out) :: Wa(Ndim,Nstat), Wb(Ndim,Nstat)

   integer :: ii, jj, kk
   real*8 :: temp

   Wa = 0.0D0; Wb = 0.0D0

   do ii=1,Nstat
   do jj=1,Sdim
      temp = Vec(jj,ii)
      ! ALPHA
      do kk=1,Ndima
         Wa(kk,ii) = Wa(kk,ii) + temp * Aa(kk,jj)
      enddo
      ! BETA
      do kk=1,Ndimb
         Wb(kk,ii) = Wb(kk,ii) + temp * Ab(kk,jj)
      enddo
   enddo
      Wa(1:Ndima,ii) = Wa(1:Ndima,ii) - Val(ii) * Ra(1:Ndima,ii)
      Wb(1:Ndimb,ii) = Wb(1:Ndimb,ii) - Val(ii) * Rb(1:Ndimb,ii)
   enddo
end subroutine open_residual
