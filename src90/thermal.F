c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              THERMAL                      \\\\\\\\\c
c=======================================================================
c
      subroutine thermal
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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit none
c
      integer  :: i, j, k, l, n, iseed
      integer :: l1, n1 
      real (rl) :: aklm, delrho(in,jn,kn), kx, ky, rmax
      real*8 :: amp, rvalue, rv, pcrbyp, TkeV
      real*8 :: mui, x, y, z, n_e, ptot, T
      real*8 :: kmin, wn, theta, phi, rphase, phase
      real*8 :: a3(in,jn,kn)
      real*8 :: gasdev, drho(in,jn,kn)
      character*9 :: inifile
      integer :: nk
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
#ifdef MPI_USED
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
#endif      
      mu = 1.0d0/(2.0*x+0.75*y+0.5625*z)
      mue = 2.0d0/(1.0d0+x)
      mui = 1.0/(1.0/mu-1.0/mue)
c 1 keV plasma; 0.1 particle /cc
      T = TkeV*1.1604e7
      d = mue*n_e*mp/mui
      p = d*boltz*T/(mu*mp)
      if (xcosmic) pcr = pcrbyp*p
      iseed = 1
      v1=0.0d0;v2=0.0d0;v3=0.0d0
      if (xmhd) then
        b1=1.0e-6;b2=1.0e-6;b3=0.0
      endif

      iseed = 1
      delrho = 0.0
      rmax = 40.*3.086e21
c
c instead of adding random perturbations, I should have the SAME perturbations for all runs. Useful to test c convergence etc.
c       
        do n = 4, 20
        do n1 = -n, n, 2*n
          kx = 2.0*pi*n1/(2.*rmax)
        do l = 4, 20
        do l1 = -l, l, 2*l
          ky = 2.0*pi*l1/(2.*rmax)
c
          call ran2(iseed,rvalue)
          phi = 2.0*pi*rvalue

          call ran2(iseed,rvalue)
          aklm = 0.15*amp*(0.5-rvalue)/sqrt(1.*n**2+l**2)

        do k = ks, ke
        do j = js, je
        do i = is, ie

          delrho(i,j,k) = delrho(i,j,k) + aklm*d(i,j,k)
     &    *cos( phi + kx*x1b(i) + ky*x2b(j) )

        enddo
        enddo
        enddo
c        
      enddo
      enddo
      enddo
      enddo

      d = d + delrho

c      write(inifile,"(a3,3i2.2)") 'bni',coords(1)
c     &                                      ,coords(2),coords(3)
c      open(unit=261,file=inifile,status='old',form='unformatted')
c      read(261) drho(is:ie,js:je,ks:ke)
c      close(261)
c      do k = ks, ke
c      do j = js, je
c      do i = is, ie
c        d(i,j,k) = d(i,j,k)*( 1.0 + amp*1.e2*drho(i,j,k))
c      enddo
c      enddo
c      enddo

      e = p/gamm1
      if (xcosmic) ecr = pcr/gamm1_cr
c
      return
      end
