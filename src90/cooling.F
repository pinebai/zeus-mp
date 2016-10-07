c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  COOLING                  \\\\\\\\\c
c                            Developed by
c                            Prateek Sharma 
c
c=======================================================================
c
       subroutine cooling 
c
c Implements f-f (+metals) cooling; assuming that CGS units are used
c Density cannot simply be scaled out like in case when there is no cooling
c or self gravity. Tozzi & Norman 2001 cooling function; see Chandran 2007, 
c Guo & Peng Oh 2007
c
c......................................................................
      use real_prec
      use config
      use param
      use root
      use grid
      use field
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use metallicity 
      use cons
      use bndry
c
      implicit NONE
c
      integer  :: i, j, k, l, n, nsubcycle
      real(rl) :: mui, n_e, n_i, TkeV
      real(rl) :: cooling_rate(in,jn,kn), const_heat, dtcool 
      real(rl) :: dvol, vol
c
      mui = 1.0d0/(1.0d0/mu - 1.0d0/mue)

      dtcool=huge

      do k=ks, ke
      do j=js, je
      do i=is, ie
        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)
c
      if (TkeV.gt.0.17235) then
          cooling_rate(i,j,k) = 1.0d-22*n_i*n_e*( 8.6d-3*TkeV**(-1.7)
     &    + 5.8d-2*TkeV**0.5 + 6.3d-2)
      else
          cooling_rate(i,j,k) = n_i*n_e*2.579e-23*(TkeV/0.17235)**6.0
      endif
c
        if (TkeV.le.1.0d-4) cooling_rate(i,j,k) = 0.0d0

        dtcool = min( dtcool, e(i,j,k)/cooling_rate(i,j,k) )
c
      enddo
      enddo
      enddo
c
#ifdef MPI_USED 
      buf_in(1) = dtcool
      call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
      dtcool  =   buf_out(1)
#endif
c
c factor of 4.0 chosen arbitrarily
c
      nsubcycle = int(4.0*dt/dtcool) + 1
c      nsubcycle = int(10.0*dt/dtcool) + 1
      dtcool = dt/nsubcycle
c
      if (mod(nhy,100).eq.0) write(17,*) nsubcycle

      do n = 1, nsubcycle
c
c keeping serial version in mind
c
      const_heat = 0.0d0
      vol = 0.0d0 
      do k = ks, ke
      do j = js, je
      do i = is, ie
        dvol = dvl1a(i)*dvl2a(j)*dvl3a(k)
        vol = vol + dvol

        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)
c
      if (TkeV.gt.0.17235) then
          cooling_rate(i,j,k) = 1.0d-22*n_i*n_e*( 8.6d-3*TkeV**(-1.7)
     &    + 5.8d-2*TkeV**0.5 + 6.3d-2)
      else
          cooling_rate(i,j,k) = n_i*n_e*2.579e-23*(TkeV/0.17235)**6.0
      endif
c
        if (TkeV.le.1.0d-4) cooling_rate(i,j,k) = 0.0d0


        const_heat = const_heat + dvol*cooling_rate(i,j,k) 
      enddo
      enddo
      enddo
#ifdef MPI_USED
      buf_in(1) = const_heat
      buf_in(2) = vol
      call MPI_ALLREDUCE( buf_in, buf_out, 2
     &                   , MPI_DOUBLE_PRECISION
     &                   , MPI_SUM, comm3d, ierr)
      const_heat = buf_out(1)
      vol = buf_out(2)
#endif
      const_heat = const_heat/vol
c
      do k = ks, ke
      do j = js, je
      do i = is, ie
c
        cooling_rate(i,j,k) = cooling_rate(i,j,k) - 0.0*const_heat

c        if (cooling_rate(i,j,k).gt.0.0) then
c          e(i,j,k) = e(i,j,k)/( 1.0 + cooling_rate(i,j,k)*dtcool/e(i,j,k) )
c        else
c          e(i,j,k) = e(i,j,k) - cooling_rate(i,j,k)*dtcool
c        endif
c
        e(i,j,k) = e(i,j,k) - dtcool*cooling_rate(i,j,k)

      enddo 
      enddo
      enddo

      enddo !subcycle
c
c making the boundary values out of date
c
      do i = 1,6
       bvstat(i,2) = 0      !  e
      enddo
c
      return
       end
