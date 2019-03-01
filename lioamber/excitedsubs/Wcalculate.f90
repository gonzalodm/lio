subroutine Wcalculate(Zvec,Dif,C,Wmat,Ndim,M)
use lrdata, only: Qvec, GxcAO, Vlr, cbas, Cocc, Cocc_trans,&
                  eigenval, root, EneSCF, Coef_trans
use garcha_mod, only: NCO
   implicit none

   integer, intent(in) :: Ndim, M
   real*8, intent(in) :: Zvec(Ndim), Dif(M,M), C(M,M)
   real*8, intent(inout) :: Wmat(M,M)

   integer :: i, j, k
   integer :: Nvirt, pos1, pos2, NCOc
   real*8 :: temp1, temp2
   real*8, dimension(:,:), allocatable :: F2e, Fxc, Ftot
   real*8, dimension(:,:), allocatable :: scratch, HXIJ, GXCIJ
   

   Nvirt = M - NCO

!  print*, "Wcalc"

!  print*, "z vector en MO basis"
!  do i=1,Ndim
!     print*, i,Zvec(i)
!  enddo
!  print*, "Gxc en AO basis"
!  do i=1,M
!  do j=1,M
!     print*, i,j,GxcAO(i,j)
!  enddo
!  enddo
!  print*,"Q vec in MO basis IA"
!  do i=1,Ndim
!     print*, i,Qvec(i)
!  enddo
!  print*,"Vec LR 1/2 en MO basis"
!  do i=1,Ndim
!     print*, i,Vlr(i)
!  enddo

   
   ! Calculate 2 electron part
   allocate(F2e(M,M)); F2e = 0.0D0
   call g2g_calculate2e(Dif,cbas,1,F2e,0)
   F2e = 2.0D0 * F2e
   print*, "F2e en AO basis"
!  do i=1,M
!  do j=1,M
!     print*, i,j,F2e(i,j)
!  enddo
!  enddo
!  stop
   
   ! Calculate xc part
   allocate(Fxc(M,M)); Fxc = 0.0D0
   call g2g_calculateg(Dif,Fxc,2)
!  print*, "Fxc en AO basis"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Fxc(i,j)
!  enddo
!  enddo
!  stop
   
   ! Obtain total 
   allocate(Ftot(M,M))
   Ftot = F2e + 2.0D0 * Fxc
   deallocate(F2e,Fxc)
!  print*, "Ftot en AO basis"
!  do i=1,M
!  do j=1,M
!     print*, i,j,Ftot(i,j)
!  enddo
!  enddo
!  stop

   allocate(scratch(M,NCO),HXIJ(NCO,NCO))
   call multlr(Ftot,Cocc,scratch,M,M,NCO,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,HXIJ,NCO,M,NCO,1.0D0,0.0D0)
   deallocate(scratch)

!  print*, "total in IJ"
!  do i=1,NCO
!  do j=1,NCO
!     print*, i,j,HXIJ(i,j)
!  enddo
!  enddo
!  stop

   ! GXC AO -> MO
   allocate(scratch(M,NCO),GXCIJ(NCO,NCO))
   call multlr(GxcAO,Cocc,scratch,M,M,NCO,1.0D0,0.0D0)
   call multlr(Cocc_trans,scratch,GXCIJ,NCO,M,NCO,1.0D0,0.0D0)
   deallocate(scratch)
!  print*, "GXC in IJ"
!  do i=1,NCO
!  do j=1,NCO
!     print*, i,j,GXCIJ(i,j)
!  enddo
!  enddo
!  stop

!  print*, "eigenvalues LR"
!  do i=1,root
!     print*, i,eigenval(i)
!  enddo
!  print*, "eigenvalues SCF"
   do i=1,M
      print*, i,EneSCF(i)
   enddo

! FORM Energy weigth in MO basis
   ! FOR BLOCK OCC x OCC
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,i
      temp1 = 0.0D0
      temp2 = 0.0D0
      print*,"guarda",NCOc-i,NCOc-j,i,j
      do k=1,Nvirt
         pos1 = (i-1)*Nvirt+k !ia
         pos2 = (j-1)*Nvirt+k !ja
         temp1 = temp1 + Vlr(pos1)*Vlr(pos2)*2.0D0
         temp2 = temp2 + EneSCF(NCO+k)*(Vlr(pos1)*Vlr(pos2)*2.0D0)
         print*, "NCO+k,pos1,pos2",NCO+k,pos1,pos2
      enddo
      !Wmat(i,j) = eigenval(root) * temp1 - temp2 &
      !            + HXIJ(i,j) + 2.0D0 * GXCIJ(i,j) 
      !Wmat(j,i) = Wmat(i,j)
      Wmat(NCOc-i,NCOc-j) = eigenval(root) * temp1 - temp2 &
                  + HXIJ(NCOc-i,NCOc-j) + 2.0D0 * GXCIJ(NCOc-i,NCOc-j) 
      Wmat(NCOc-j,NCOc-i) = Wmat(NCOc-i,NCOc-j)
      !print*, "IJ",i,j,Wmat(NCOc-i,NCOc-j)
    enddo
    enddo
 
    ! FORM BLOCK VIR x VIR
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
          print*,"NCOc-k,pos1,pos2",NCOc-k,pos1,pos2
          !print*, "k",EneSCF(NCOc-k),Vlr(pos1),Vlr(pos2)
       enddo
       Wmat(NCO+i,NCO+j) = eigenval(root)*temp1 + temp2
       Wmat(NCO+j,NCO+i) = Wmat(NCO+i,NCO+j)
       print*, "AB",NCO+i,NCO+j,Wmat(NCO+i,NCO+j)
    enddo
    enddo
       
    ! FORM BLOCK OCC x VIR
    NCOc = NCO + 1
    do i=1,NCO
    do j=1,Nvirt
       pos1 = (i-1)*Nvirt+j !ij
       Wmat(NCO+j,NCOc-i) = Qvec(pos1) + EneSCF(NCOc-i)*Zvec(pos1)
       Wmat(NCOc-i,NCO+j) = Wmat(NCO+j,NCOc-i)
       print*,"NCOc-i,pos1",NCOc-i,pos1
       print*, "IA",NCO+j,NCOc-i,Wmat(NCO+j,NCOc-i)
    enddo
    enddo
    !stop

    !do i=1,M
    !   Wmat(i,i) = Wmat(i,i) * 0.5D0
    !enddo

    print*, "W in MO"
    do i=1,M
    do j=1,M
       print*, i,j,Wmat(i,j)
    enddo
    enddo
    ! Wmat MO -> AO
    allocate(scratch(M,M))
    call multlr(C,Wmat,scratch,M,M,M,1.0D0,0.0D0)
    call multlr(scratch,Coef_trans,Wmat,M,M,M,1.0D0,0.0D0)
    deallocate(scratch)
    print*, "W in AO"
    do i=1,M
    do j=1,M
       print*, i,j,Wmat(i,j)
    enddo
    enddo
end subroutine Wcalculate
