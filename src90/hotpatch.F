c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              HOTPATCH                     \\\\\\\\\c
c=======================================================================
c
      subroutine hotpatch
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
      use root
      use domain
      use mpiyes
      use mpipar
c
      implicit none
c
      integer  :: i, j, k
      real*8 :: pcrbyp, TkeV, amp
      real*8 :: mui, x, y, z, n_e, ptot, T
      namelist  /pgen/ x, y, z, amp, pcrbyp, n_e, TkeV
c
c-----------------------------------------------------------------------
c
      n_e = 0.1
      TkeV = 1.0
      pcrbyp = 0.001
      x = 0.7
      y = 0.28
      z = 0.02
      amp = 0.01
      if (myid.eq.0) then
        read(1,pgen)
        write(2,pgen)
        buf_in(1) = x
        buf_in(2) = y
        buf_in(3) = z
        buf_in(4) = amp
        buf_in(5) = pcrbyp
        buf_in(6) = n_e
        buf_in(7) = TkeV
      endif
        call MPI_BCAST( buf_in, 7, MPI_DOUBLE_PRECISION
     &                , 0, comm3d, ierr )
      if (myid.ne.0) then
        x = buf_in(1)
        y = buf_in(2)
        z = buf_in(3)
        amp = buf_in(4)
        pcrbyp = buf_in(5)
        n_e = buf_in(6)
        TkeV = buf_in(7)
      endif
      mu = 1.0d0/(2.0*x+0.75*y+0.5625*z)
      mue = 2.0d0/(1.0d0+x)
      mui = 1.0/(1.0/mu-1.0/mue)
c 1 keV plasma; 0.1 particle /cc
      T = TkeV*1.1604e7
      d = mue*n_e*mp/mui
      p = d*boltz*T/(mu*mp)
      if (xcosmic) pcr = pcrbyp*p
      v1=0.0d0;v2=0.0d0;v3=0.0d0
      if (xmhd) then
        b1=1.0e-6;b2=1.0e-6;b3=0.0
      endif
      do k = 1, kn
      do j = 1, jn 
      do i = 1, in
        if (x1b(i).ge.10.*3.086e21.and.x1b(i).le.20.0*3.086e21
     &  .and.x2b(j).ge.10.*3.086e21.and.x2b(j).le.20.0*3.086e21) 
     &  d(i,j,k) = 0.5*d(i,j,k)
c        if (x1b(i).ge.10.*3.086e21.and.x1b(i).le.20.0*3.086e21
c     &  .and.x2b(j).ge.20.*3.086e21.and.x2b(j).le.30.0*3.086e21)
c     &  d(i,j,k) = 0.5*d(i,j,k)
c        if (x1b(i).ge.20.*3.086e21.and.x1b(i).le.30.0*3.086e21
c     &  .and.x2b(j).ge.20.*3.086e21.and.x2b(j).le.30.0*3.086e21)
c     &  d(i,j,k) = 0.5*d(i,j,k)
c        if (x2b(j).ge.10.*3.086e21.and.x2b(j).le.20.0*3.086e21)
c     &  d(i,j,k) = 0.5*d(i,j,k)
c        if (x1b(i).ge.10.*3.086e21.and.x1b(i).le.20.0*3.086e21)
c     &  d(i,j,k) = 0.5*d(i,j,k)
      enddo
      enddo
      enddo
      e = p/gamm1
      if (xcosmic) ecr = pcr/gamm1_cr
c
      return
      end
