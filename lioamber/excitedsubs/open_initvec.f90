subroutine open_vec_init(vecA,vecB,Ea,Eb,M,NCOa,NCOb,N,vecnum)
   implicit none

   integer, intent(in) :: M, NCOa, NCOb, N, vecnum
   real*8, intent(in) :: Ea(M), Eb(M)
   real*8, intent(out) :: vecA(N,vecnum), vecB(N,vecnum)

   integer :: i, num, Ndima, Ndimb, ind, occ, virt, NCOc, Nvirt
   integer, dimension(:), allocatable :: indexs
   real*8, dimension(:), allocatable :: deltaE

   vecA = 0.0D0
   vecB = 0.0D0
  
   Ndima = NCOa * (M - NCOa)
   Ndimb = NCOb * (M - NCOb)
   allocate(indexs(Ndima+Ndimb),deltaE(Ndima+Ndimb))

!  CALCULATE ENERGIES DIFFERENCE
   ! alpha
   NCOc = NCOa + 1
   Nvirt = M - NCOa
   do i=1,Ndima
     ind = i - 1
     occ = NCOa - (ind/Nvirt)
     virt = mod(ind,Nvirt) + NCOc
     deltaE(i) = Ea(virt) - Ea(occ)
   enddo

   ! beta
   NCOc = NCOb + 1
   Nvirt = M - NCOb
   do i=1,Ndimb
     ind = i - 1
     occ = NCOb - (ind/Nvirt)
     virt = mod(ind,Nvirt) + NCOc
     deltaE(Ndima+i) = Eb(virt) - Eb(occ)
   enddo

!  SORTING OF ENERGIES DIFERENCES
   call iargsort(deltaE,indexs,Ndima+Ndimb)

!  INITIAL TRIALS VECTORS
   do i=1,vecnum
      if (indexs(i) <= Ndima) then
         vecA(indexs(i),i) = 1.0D0
      else
         vecB(indexs(i)-Ndima,i) = 1.0D0
      endif
   enddo
   deallocate(indexs,deltaE)

end subroutine open_vec_init

subroutine iargsort(a,b,N)

   implicit none

   integer, intent(in) :: N
   real*8, intent(in):: a(N)
   integer, intent(out) :: b(size(a))
   integer :: i,imin
   real*8 :: temp
   real*8 :: a2(size(a))

   a2 = a
   do i = 1, N
      b(i) = i
   end do
   do i = 1, N-1
      imin = minloc(a2(i:),1) + i - 1
      if (imin /= i) then
         temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
         temp = b(i); b(i) = b(imin); b(imin) = temp
      end if
   end do
end subroutine
