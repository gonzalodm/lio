subroutine vecMOtomatMO(vec,mat,M,NCO,Nvirt,Sdim,Nvec,V1,dim)
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Nvec, V1, dim, Sdim
   real*8, intent(in) :: vec(dim,Sdim)
   real*8, intent(out) :: mat(M,M,Nvec)

   integer :: i,row,col,NCOc

   NCOc = NCO - 1
   mat = 0.0D0

   do i=1,Nvec
     do row=0,NCOc
     do col=1,Nvirt
        mat(NCOc-row+1,NCO+col,i) = vec(row*Nvirt+col,V1+i)
     enddo
     enddo
   enddo
end subroutine vecMOtomatMO

