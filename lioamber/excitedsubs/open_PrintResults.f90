subroutine open_PrintResults(VecA,VecB,Val,O,Ndim,Nstat,M,NCOa,NCOb)
use lrdata, only: nfo
   implicit none

   integer, intent(in) :: Ndim, Nstat, M, NCOa, NCOb
   real*8, intent(in) :: VecA(Ndim,Nstat), VecB(Ndim,Nstat), Val(Nstat), O(Nstat)

   character(len=4) :: j_char, from_char, to_char
   integer :: ii, jj, from, to, Ndima, Ndimb
   real*8 :: value_x

   Ndima = (M - NCOa) * NCOa
   Ndimb = (M - NCOb) * NCOb
   nfo = 0

   do jj=1,Nstat
     write (j_char, '(i4)') jj
     write(*,100) adjustl(j_char), Val(jj), 45.56335D0/val(jj), O(jj)
     ! ALPHA
     from = NCOa
     to = NCOa + 1
     do ii=1,Ndima
        value_X = VecA(ii,jj)
        if ( abs(value_X) > 0.1D0 ) then
         write (from_char, '(i4)') from+nfo
         write (to_char, '(i4)') to+nfo
         write(*,101) adjustl(from_char), adjustl(to_char), "(A)", value_X
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
        value_X = VecB(ii,jj)
        if ( abs(value_X) > 0.1D0 ) then
         write (from_char, '(i4)') from+nfo
         write (to_char, '(i4)') to+nfo
         write(*,101) adjustl(from_char), adjustl(to_char), "(B)", value_X
        endif
        to = to + 1
        if ( to == M+1 ) then
            from = from - 1
            to = NCOb + 1
        endif
     enddo
     print*, " "
   enddo
   
   100 FORMAT(1X,"STATE ",A,3X,"ENERGY=",F8.4," Hartree, ",&
              F12.6," nm"," OSC=",F8.4)
   101 FORMAT(6X,A," -> ",A,1X,A,2X,F14.7)
end subroutine open_PrintResults
