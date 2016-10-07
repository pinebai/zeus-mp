c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 RTI                       \\\\\\\\\c
c=======================================================================
c
      subroutine rti
c
c Rayleigh Taylor instability in cartesian geometry with gravity in $-x1$ direction
c acceleration due to grav. is taken to be 1.0
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use grid
      use cons
      use metallicity
      use gravmod
      use root
      use domain
      use mpiyes
      use mpipar
c
      implicit none
c
      integer  :: i, j, k, l, iseed
      real*8 :: ptot, amp, rvalue
      namelist  /pgen/ amp
c
c-----------------------------------------------------------------------
c
      ptot = 10.0
      amp = 0.001
      iseed = 1
      do k = 1, kn
      do j = 1, jn
      do i = 1, in
        call ran2(iseed,rvalue)
        v1(i,j,k) = amp*(rvalue-0.5)*sqrt(ptot)
        call ran2(iseed,rvalue)
        v2(i,j,k) = amp*(rvalue-0.5)*sqrt(ptot)
        call ran2(iseed,rvalue)
        v3(i,j,k) = amp*(rvalue-0.5)*sqrt(ptot)
        if (xmhd) then
          b1(i,j,k) = 0.0d0
          b2(i,j,k) = 1.0d0
          b3(i,j,k) = 0.0d0
        endif
        if (x1b(i).lt.0.0) then
          d(i,j,k) = 0.1
        else
          d(i,j,k) = 1.0
        endif 
      enddo
      enddo
      enddo
c
      p(ie+3,:,:) = ptot
c
      do l = ntiles(1)-1, 0, -1
        if (l.eq.coords(1)) then
          if (l.ne.ntiles(1)-1) then
            nreq=nreq+1
            call MPI_RECV(ptot, 1, MPI_DOUBLE_PRECISION, n1p, 12345
     &      , comm3d, stat, ierr)
            call MPI_WAIT(stat, req(nreq), ierr)   
            p(ie+3,:,:) = ptot 
          endif
          do k = 1, kn
          do j = 1, jn
          do i = in-1, 1, -1
c
c g chosen to be 1.0
c
            p(i,j,k) = p(i+1,j,k) 
     &               + 0.5*(d(i,j,k)+d(i-1,j,k))*1.0*dx1b(i+1)
          enddo
          enddo
          enddo     
          ptot = p(is+2,(js+je)/2,(ks+ke)/2) 
          if (l.ne.0) then
            nreq = nreq + 1
            call MPI_ISEND(ptot, 1, MPI_DOUBLE_PRECISION, n1m, 12345
     &      , comm3d, req(nreq), ierr)
          endif
        endif
      enddo 
c
      e = p/gamm1
c
      if (xcosmic) then
        j=(js+je)/2
        k=(ks+ke)/2 
        do i=1, in
          if (x1b(i).lt.0.0) then
            pcr(i,:,:) = p(i,:,:)*0.99
           else
            pcr(i,:,:) = p(i,:,:)*0.01
          endif
          p(i,:,:) = p(i,:,:) - pcr(i,:,:)
          ecr(i,:,:) = pcr(i,:,:)/(gamma_cr-1.0)
          e(i,:,:) = p(i,:,:)/(gamma-1.0)
        enddo
      endif
c
c      j=(js+je)/2
c      k=(ks+ke)/2
c      do i=1, in
c        write(22,1234) x1b(i), d(i,j,k), p(i,j,k), pcr(i,j,k)
c      enddo
c
c1234  format(4e20.7)
c
c      stop
c
      return
      end
c====================================================================
c  numerical recipes random number generator ran2
c    requires input seed value=iseed
c    returns real random number=rvalue
c    Also updates iseed for next call 
c
      subroutine ran2(iseed,rvalue)
      integer iseed
      real*8 rvalue 
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,jj,kk,iv(NTAB),iy
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
c
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
        iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum 
      return
      end
