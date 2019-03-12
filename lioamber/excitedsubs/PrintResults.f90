subroutine PrintResults(vec,val,O,N,nstat,Mlr,NCOlr)
use lrdata, only: nfo
   implicit none

   integer, intent(in) :: N, nstat, Mlr, NCOlr
   real*8, intent(in) :: vec(N,nstat),val(nstat),O(nstat)

   integer :: i,j,from,to
   real*8 :: value_X

   from = NCOlr
   to = NCOlr + 1

   do j=1, nstat
   write(*,100) j, val(j), 45.56335D0/val(j), O(j)
   do i=1, N
      value_X = vec(i,j) / dsqrt(2.0D0)
      if ( abs(value_X) > 0.1D0 ) then
         write(*,101) from+nfo, to+nfo, value_X
      endif
      to = to + 1
      if ( to == Mlr+1 ) then
          from = from - 1
          to = NCOlr + 1
      endif
   enddo
      print*, " "
      from = NCOlr
      to = NCOlr + 1
   enddo

   100 FORMAT(1X,"STATE ",I2,3X,"ENERGY=",F8.4," Hartree, ",&
              F12.6," nm"," OSC=",F8.4)
   101 FORMAT(6X,I2," -> ",I2,2X,F14.7)
end subroutine PrintResults
