subroutine fcaApp(C,E,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use lrdata, only: FCA, nfo, nfv
   implicit none
   
   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   real*8, allocatable, intent(inout) :: C(:,:), E(:)

   integer :: i
   real*8, dimension(:,:), allocatable :: C_orig
   real*8, dimension(:), allocatable :: E_orig

   if (FCA) then
      print*, "Using Frozen Core and Valence Approximation"

      allocate(C_orig(M,M),E_orig(M))
      C_orig = C; E_orig = E
      
      deallocate(C,E)
      Nvirt = M - NCO - nfv
      NCOlr = NCO - nfo
      Mlr = M - nco - nfv
      Ndim = Nvirt * NCOlr
      allocate(C(M,Mlr),E(Mlr))

      do i=1, NCOlr
        C(:,i) = C_orig(:,i+nfo)
        E(i) = E_orig(i+nfo)
      enddo
      do i=1, Nvirt
         C(:,NCOlr+i) = C_orig(:,i+NCO)
         E(NCOlr+i) = E_orig(i+NCO)
      enddo
      deallocate(C_orig,E_orig)
   else
      print*, "FCA is not used"
      nfo = 0
      Nvirt = M - NCO
      NCOlr = NCO
      Mlr = M
      Ndim = Nvirt * NCOlr
   endif
end subroutine fcaApp
