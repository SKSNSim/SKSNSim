**********************************************************************
*     ---------------------------------------- 
      subroutine read_timevent(runnum, totevents)
*     ----------------------------------------
*
*     ( purpose )
*      Calculate the number of event from timevent file( SRN NEUTRINO )
*
*     ( input )
*       nrun    : run date
*
*     ( output )
*       totevents : number of event proportional to live time
*
*     ( creation date and author )
*       2022.05.24 ;  M.Harada
*
**********************************************************************
      implicit none


      INTEGER IERR,runnum, IDATE(3),ITIME(3),istat, totevents
      INTEGER iwrite,lfbadrun, i
      INTEGER ibad_subrun(3000)
      INTEGER idate_prev(3)/0, 0, 0/
      INTEGER itime_prev(3)/0, 0, 0/ 
      REAL tweight, sdir(3), rnd, fracleft
      REAL totaltime
      character*200 timefile, fname
      INTEGER nrun,nsub
      logical*4 ex
      integer lfirst
      integer*4 nevents, nev, ievent
      integer num
      parameter(num = 1) ! Number of events per minute
      real rlu

      integer lenchr

      if ( (runnum .ge. 60000) .and.  (runnum .lt. 79999) ) then 
        timefile =
     $    '/home/sklowe/realtime_sk4_rep/solar_apr19/timevent/timevent.r061525.r077958'
      else if ( (runnum .ge. 80000) .and. (runnum .lt. 84999) ) then  
        timefile = 
     $    '/home/sklowe/realtime_sk5_rep/solar_nov20/timevent/timevent.r080539.r082915'
      else if ( runnum .ge. 85000 ) then 
      print *, runnum
        timefile = 
     $    '/home/sklowe/realtime_sk6_rep/solar_may22/timevent/timevent.r080539.r087073'
        print *, timefile
      endif

      CALL KZINIT

      lfirst = 0
c*** open time event file
      print*,'Open timefile:',timefile
      open(10, file=timefile,
     $     status='old',form='unformatted')
      print*,'timefile opened'

      print*,'Start Read '
 10   continue
****  read run number, subrun number, date(yy,mm,dd), time(hh,mm,ss)
****  run time(sec), solar direction(direction cosine toward(?) sun), istat
****  '
      call slredtimev(10,nrun,nsub,idate,itime,tweight,sdir,istat)
      nevents = 0
 
      if( nrun .lt. 60000 )then
        nrun = nrun + 65536         ! for SK-IV analysis
      endif
c***      print *, nrun
 
      if(istat.ne.0) then !istat=-1 EOF 
          if(lfirst .eq. 1) goto 100  ! at least one sub run is read
          print *, 'reading error! stop'
          stop
      endif

      if(nrun.ne.runnum) goto 10

      lfirst = 1

c      if( nrun .lt. 60000 )then
c          nrun = nrun + 65536         ! for SK-IV analysis
c      endif


c*** check duplicate date
      if (idate(1).eq.idate_prev(1) .and.
     $     idate(2).eq.idate_prev(2) .and.
     $     idate(3).eq.idate_prev(3) .and.
     $     itime(1).eq.itime_prev(1) .and.
     $     itime(2).eq.itime_prev(2) .and.
     $     itime(3).eq.itime_prev(3)) then

          if (nsub.gt.0 .and. nsub.lt.3000) then
              ibad_subrun(nsub) = 1
          else
              print *, 'data error!', nrun, nsub,idate,itime
              stop
          endif
      endif

c*** save previous date
      idate_prev(1) = idate(1)
      idate_prev(2) = idate(2)
      idate_prev(3) = idate(3)
      itime_prev(1) = itime(1)
      itime_prev(2) = itime(2)
      itime_prev(3) = itime(3)


c***  badsub run cut
      if( lfbadrun(nrun, nsub).eq.1) then
c            print *, "Bad run!"
          goto 10
      endif

c** multiple badsub run cut
      if (nsub.gt.0 .and. nsub.lt.3000) then
          if (ibad_subrun(nsub).ne.0) then
c              print *, "Bad subrun!"
              goto 10
          endif
      endif

      ievent = ievent + 1
c      if(mod(ievent,10000).eq.0) print*,ievent,' is finish'

c*** make vector ("num" times per minute)
      nev = int(tweight/60. * num) + 1
      fracleft = tweight/60. * num - nev + 1

cccccccccccccccccccccccccccccccccccccccccccccccc
      DO I=1,nev
c       If number of events is fractional, complete by drawing random number
        if (I.eq.nev) then
          rnd = rlu()
          if (rnd.gt.fracleft) exit
        endif
        nevents = nevents + 1
      enddo

      totaltime = totaltime + tweight
      totevents = totevents + nevents
cccccccccccccccccccccccccccccccccccccccccccccccc
c      print *, "total ", tweight, nevents, nev, totaltime, totevents
      goto 10

100   continue
      close(10)
      print *, "total :", totaltime, totevents

      return
      END
