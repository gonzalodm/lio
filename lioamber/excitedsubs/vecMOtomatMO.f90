subroutine vecMOtomatMO(vec,mat,M,NCO,Nvirt,Sdim,Nvec,V1,dim)
use lrdata, only: ppTDA
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Nvec, V1, dim, Sdim
   real*8, intent(in) :: vec(dim,Sdim)
   real*8, intent(out) :: mat(M,M,Nvec)

   integer :: i,row,col,NCOc
   integer :: count_ppTDA

   NCOc = NCO - 1
   mat = 0.0D0

   if (.not. ppTDA) then ! TDA

     do i=1,Nvec
       do row=0,NCOc
       do col=1,Nvirt
          mat(NCOc-row+1,NCO+col,i) = vec(row*Nvirt+col,V1+i)
       enddo
       enddo
     enddo

   else ! pp TDA

     do i=1,Nvec
       count_ppTDA = 1
       do row=1,Nvirt
       do col=row+1,Nvirt
          mat(NCO+row,NCO+col,i) = vec(count_ppTDA,V1+i)
          !print*, NCO+row,NCO+col,count_ppTDA, vec(count_ppTDA,V1+i)
          count_ppTDA = count_ppTDA + 1
       enddo
       enddo
       print*, ""
     enddo

   endif
end subroutine vecMOtomatMO

