subroutine Wcalculate(Zvec,Dif,Qvec,GxcAO,Vlr,C,dE,EneSCF,Wmat,Ndim,M,NCO)
use lrdata, only: cbas, Cocc, Cocc_trans, Coef_trans, fitLR
   implicit none

   integer, intent(in) :: Ndim, M, NCO
   real*8, intent(in) :: Zvec(Ndim), Dif(M,M), Qvec(Ndim)
   real*8, intent(in) :: GxcAO(M,M), Vlr(Ndim), C(M,M), EneSCF(M)
   real*8, intent(in) :: dE
   real*8, intent(inout) :: Wmat(M,M)

   integer :: i, j, k
   integer :: Nvirt, pos1, pos2, NCOc
   real*8 :: temp1, temp2
   real*8, dimension(:,:), allocatable :: F2e, Fxc, Ftot, Dcopy
   real*8, dimension(:,:), allocatable :: scratch, HXIJ, GXCIJ
   
   Nvirt = M - NCO
   ! Calculate 2 electron part
   allocate(F2e(M,M)); F2e = 0.0D0

   if (.not. fitLR) then
      call g2g_calculate2e(Dif,cbas,1,F2e,1)
      F2e = 2.0D0 * F2e
   else 
      allocate(Dcopy(M,M)); Dcopy = Dif
      call calc2eFITT(Dcopy,F2e,1,M)
      deallocate(Dcopy)
   endif
   
   ! Calculate xc part
   allocate(Fxc(M,M)); Fxc = 0.0D0
   call g2g_calculateg(Dif,Fxc,2)
   
   ! Obtain total 
   allocate(Ftot(M,M))
   Ftot = F2e + 2.0D0 * Fxc
   deallocate(F2e,Fxc)

   ! Change basis Fock Total. AO - > MO(OCC x OCC)
   allocate(scratch(M,NCO),HXIJ(NCO,NCO))
   call multlr(Ftot,Cocc,scratch,M,M,NCO,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,HXIJ,NCO,M,NCO,1.0D0,0.0D0)
   deallocate(scratch,Ftot)

   ! Change basis GXC. AO -> MO(OCC x OCC)
   allocate(scratch(M,NCO),GXCIJ(NCO,NCO))
   call multlr(GxcAO,Cocc,scratch,M,M,NCO,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,GXCIJ,NCO,M,NCO,1.0D0,0.0D0)
   deallocate(scratch)

!  FORM ENERGY WEIGTHED DIFFERENCE DENSITY MATRIX
   ! BLOCK OCC x OCC
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,i
      temp1 = 0.0D0
      temp2 = 0.0D0
      do k=1,Nvirt
         pos1 = (i-1)*Nvirt+k !ia
         pos2 = (j-1)*Nvirt+k !ja
         temp1 = temp1 + Vlr(pos1)*Vlr(pos2)*2.0D0
         temp2 = temp2 + EneSCF(NCO+k)*(Vlr(pos1)*Vlr(pos2)*2.0D0)
      enddo
      Wmat(NCOc-i,NCOc-j) = dE * temp1 - temp2 &
                  + HXIJ(NCOc-i,NCOc-j) + 2.0D0 * GXCIJ(NCOc-i,NCOc-j) 
      Wmat(NCOc-j,NCOc-i) = Wmat(NCOc-i,NCOc-j)
   enddo
   enddo
   deallocate(HXIJ,GXCIJ)
 
   ! BLOCK VIR x VIR
   NCOc = NCO + 1
   do i=1,Nvirt
   do j=1,i
      temp1 = 0.0D0
      temp2 = 0.0D0
      do k=1,NCO
         pos1 = (k-1)*Nvirt+i !ki
         pos2 = (k-1)*Nvirt+j !kj
         temp1 = temp1 + Vlr(pos1)*Vlr(pos2)*2.0D0
         temp2 = temp2 + EneSCF(NCOc-k)*(Vlr(pos1)*Vlr(pos2)*2.0D0)
      enddo
      Wmat(NCO+i,NCO+j) = dE * temp1 + temp2
      Wmat(NCO+j,NCO+i) = Wmat(NCO+i,NCO+j)
   enddo
   enddo
       
   ! BLOCK OCC x VIR
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,Nvirt
      pos1 = (i-1)*Nvirt+j !ij
      Wmat(NCO+j,NCOc-i) = Qvec(pos1) + EneSCF(NCOc-i)*Zvec(pos1)
      Wmat(NCOc-i,NCO+j) = Wmat(NCO+j,NCOc-i)
   enddo
   enddo

   ! Change basis Wmat. MO -> AO
   allocate(scratch(M,M))
   call multlr(C,Wmat,scratch,M,M,M,1.0D0,0.0D0)
   call multlr(scratch,Coef_trans,Wmat,M,M,M,1.0D0,0.0D0)
   deallocate(scratch)
end subroutine Wcalculate
