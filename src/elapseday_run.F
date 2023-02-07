!===============================================
      subroutine elapseday_run(r,ed)
!===============================================
        implicit none
#include "skday.h"
        integer r,ed
        integer i_run
        real ed_r

        call skrunday


        ed = int(relapse(r))

      end subroutine elapseday_run
