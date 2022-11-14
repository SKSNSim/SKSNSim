! interface to skelapseday.F defined in SKOFL
      subroutine elapseday_date(y,m,d,ed)
        implicit none
#include "skday.h"
        integer y,m,d,ed
        integer i_run
        real ed_r

        call skrunday

        ed_r = 999999 ! skelapseday(nday,ntim)
        do i_run = 65000,999999
!           print *,i_run
!           print *,i_run,ryear(i_run),rmon(i_run),rday(i_run),relapse(i_run)
           if( (ryear(i_run).eq.y) .and.
     $         (rmon(i_run).eq.m) .and.
     $         (rday(i_run).le.d .and. d.le.rday(i_run))) then
               ed_r = relapse(i_run)
               exit
           end if
        end do
        ed = int(ed_r)

      end subroutine elapseday_date
