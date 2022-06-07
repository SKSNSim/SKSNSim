      subroutine slredtimev(ilun,nrun,nsub,id,it,t,sdir,istat)
c
c   read time event
c
c     input :   ilun     logical unit number (-1: reset idata, ndata)
c     output:   nrun     run #
c               nsub     subrun #
c               id(3)    date
c               it(3)    time
c               t        weight of the event (sec)
c               sdir(3)  direction pointing from the sun
c               istat    status code  =0 OK   =-1 end of file
c
c     10-NOV-2000 Y.Takeuchi   support ilun = -1 (reset idata, ndata)
c     22-OCT-2007 Y.Takeuchi   ishft,iand are used for linux
c
      dimension id(3),it(3),sdir(3)
      dimension nruns(1000),istime(1000),ts(1000),s(3,1000)
      data ndata/0/
      data idata/0/
      save idata,ndata
      istat=0
ccc add by YT cccccccccccccc
      if (ilun.lt.0) then
         ndata = 0
         idata = 0
         print *, "slredtimev: initialized"
         return
      endif
cccccccccccccccccccccccccccc
      if(ndata.eq.0 .or. idata.eq.ndata)then
        read(ilun,end=1000)ndata,(nruns(i),istime(i),ts(i),
     &                            s(1,i),s(2,i),s(3,i),i=1,ndata)
c        do i=0, ndata
c         print *, "****test****", ndata, ishft(nruns(i), -16),
c     &   iand(nruns(i), 65535)
c         enddo
        idata=0
      endif
      idata=idata+1
c      nrun=nruns(idata)/65536
c      nsub=mod(nruns(idata),65536)
      nrun=ishft(nruns(idata),-16)
      nsub=iand(nruns(idata),65535)
      call slidit(id,it,istime(idata))
      t=ts(idata)
      sdir(1)=s(1,idata)
      sdir(2)=s(2,idata)
      sdir(3)=s(3,idata)
      return
 1000 istat=-1
      return
      end
