subroutine open_PrintResults(VecA,VecB,Val,O,Ndim,Nstat,M,NCOa,NCOb)
use lrdata, only: nfo
   implicit none

   integer, intent(in) :: Ndim, Nstat, M, NCOa, NCOb
   real*8, intent(in) :: VecA(Ndim,Nstat), VecB(Ndim,Nstat), Val(Nstat), O(Nstat)

   integer :: ii, jj, from, to, Ndima, Ndimb
   real*8 :: value_x

   Ndima = (M - NCOa) * NCOa
   Ndimb = (M - NCOb) * NCOb
   nfo = 0

   do jj=1,Nstat
     write(*,100) jj, Val(jj), 45.56335D0/Val(jj), O(jj)
     ! ALPHA
     from = NCOa
     to = NCOa + 1
     do ii=1,Ndima
        if ( abs(VecA(ii,jj)) > 0.1D0 ) then
           write(*,101) from+nfo, to+nfo, "(A)",VecA(ii,jj)
        endif
        to = to + 1
        if ( to == M+1 ) then
            from = from - 1
            to = NCOa + 1
        endif
     enddo
     ! BETA
     from = NCOb
     to = NCOb + 1
     do ii=1,Ndimb
        if ( abs(VecB(ii,jj)) > 0.1D0 ) then
           write(*,101) from+nfo, to+nfo, "(B)",VecB(ii,jj)
        endif
        to = to + 1
        if ( to == M+1 ) then
            from = from - 1
            to = NCOb + 1
        endif
     enddo
     print*, " "
   enddo
   
   100 FORMAT(1X,"STATE ",I2,3X,"ENERGY=",F8.4," Hartree, ",&
              F12.6," nm"," OSC=",F8.4)
   101 FORMAT(6X,I3," -> ",I3,1X,A,2X,F14.7)
end subroutine open_PrintResults
