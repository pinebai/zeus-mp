c=======================================================================
c
      subroutine printd 
c
c  PURPOSE: dumps history file; conserved volume integrated quantities, etc.
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer, parameter :: nscal = 30 
      integer :: i, j, k
      real(rl) :: scal(nscal), dvol, darea, dm 
      real(rl) :: dmin, dmax, Tmin
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      scal = 0.0d0 
      scal(1) = time
      scal(2) = dt
      dmin = huge
      dmax = tiny
      Tmin = huge
      do k = ks, ke
      do j = js, je
      do i = is, ie
        dvol = dvl1a(i)*dvl2a(j)*dvl3a(k)
        dm = dvol*d(i,j,k)
c total mass
        scal( 3) = scal( 3) + dm
c total internal energy
        scal( 4) = scal( 4) + dvol*e(i,j,k)
c total (r,t,p) KE
        scal( 5) = scal( 5) + dm*0.25*(v1(i,j,k)**2+v1(i+1,j,k)**2)
        scal( 6) = scal( 6) + dm*0.25*(v2(i,j,k)**2+v2(i,j+1,k)**2)
        scal( 7) = scal( 7) + dm*0.25*(v3(i,j,k)**2+v3(i,j,k+1)**2)
c total (r,t,p) magnetic energy
        if (xmhd) then
          scal( 8) = scal( 8) + dvol*0.25*(b1(i,j,k)**2+b1(i+1,j,k)**2)
          scal( 9) = scal( 9) + dvol*0.25*(b2(i,j,k)**2+b2(i,j+1,k)**2)
          scal(10) = scal(10) + dvol*0.25*(b3(i,j,k)**2+b3(i,j,k+1)**2)
        endif
c gravitational potential energy; for now, either xptmass or xdmpotl
        if (xptmass) then
          scal(11) = scal(11) - (guniv*ptmass/(x1b(i)-rg))*dm
        endif
        if (xdmpotl) then
          scal(11) = scal(11) + phi_dm(i)*dm
        endif
c minimum density
          scal(16) = min(dmin,d(i,j,k))
c maximum density
          scal(18) = max(dmax,d(i,j,k))
c minimum temperature
          scal(17) = min(Tmin,e(i,j,k)/d(i,j,k))
c angular momentum
          scal(14) = scal(14) + dm*0.5*(v3(i,j,k)+v3(i,j,k+1))
     &             *g31b(i)*g32b(j)
c cosmic ray energy
        if (xcosmic) then
          scal(15) = scal(15) + dvol*ecr(i,j,k)
        endif
      enddo
      enddo
      enddo
      scal(12) = eadd
      scal(13) = madd    
      buf_in(1:13) = scal(3:15)
      buf_in(14:15) = scal(16:17)
      buf_in(16) = scal(18)
#ifdef MPI_USED
c
c MPI_2DOUBLE_PRECISION was not working fine with MPI_SUM
c
      call MPI_REDUCE(buf_in, buf_out, 13, MPI_DOUBLE_PRECISION
     &, MPI_SUM, 0, comm3d, ierr)
      call MPI_REDUCE(buf_in(14), buf_out(14), 2, MPI_DOUBLE_PRECISION
     &, MPI_MIN, 0, comm3d, ierr)
      call MPI_REDUCE(buf_in(16), buf_out(16), 1, MPI_DOUBLE_PRECISION
     &, MPI_MAX, 0, comm3d, ierr)
      if (myid.eq.0) then
        scal(3:15) = buf_out(1:13)
        scal(16:17) = buf_out(14:15) 
        scal(18) = buf_out(16)
        write(3,2001) (scal(i), i=1, 18)
      endif
#else
        write(3,2001) (scal(i), i=1, 18)
#endif
2001  format(18e12.4)
c radial luminosities through inner boundary; integrated over theta and phi
      i = is
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(19) = scal(19) + dm
c internal energy
        scal(20) = scal(20) + darea*e(i,j,k)*v1(i,j,k)
c angular momentum
        scal(21) = scal(21) + dm*v3(i,j,k)*g31a(i)*g32b(j)
c kinetic energy
        scal(22) = scal(22) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
c Poynting E*B
        if (xmhd) then
          scal(23) = scal(23) +  darea*(
     &             v3(i,j,k)*b1(i,j,k)*b3(i,j,k)
     &             +0.25*(v2(i,j,k)+v2(i,j+1,k))*b1(i,j,k)
     &             *(b2(i,j,k)+b2(i,j+1,k))
     &             -v1(i,j,k)*(b3(i,j,k)**2
     &             +0.5*(b2(i,j,k)**2+b2(i,j+1,k)**2)))
        endif
        if (xcosmic) then
          scal(29) = scal(29) + darea*ecr(i,j,k)*v1(i,j,k)
        endif
      enddo
      enddo
c radial luminosities through outer boundary; integrated over theta and phi
      i = ie+1
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(24) = scal(24) + dm
c internal energy
        scal(25) = scal(25) + darea*e(i,j,k)*v1(i,j,k)
c angular momentum
        scal(26) = scal(26) + dm*v3(i,j,k)*g31a(i)*g32b(j)
c kinetic energy
        scal(27) = scal(27) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
c Poynting E*B
        if (xmhd) then
          scal(28) = scal(28) +  darea*(
     &             v3(i,j,k)*b1(i,j,k)*b3(i,j,k)
     &             +0.25*(v2(i,j,k)+v2(i,j+1,k))*b1(i,j,k)
     &             *(b2(i,j,k)+b2(i,j+1,k))
     &             -v1(i,j,k)*(b3(i,j,k)**2
     &             +0.5*(b2(i,j,k)**2+b2(i,j+1,k)**2)))
        endif
        if (xcosmic) then
          scal(30) = scal(30) + darea*ecr(i,j,k)*v1(i,j,k)
        endif
      enddo
      enddo
      write(51,2002) time, (scal(i), i=19, 30)
2002  format(13e12.4)
c
      return
      end
