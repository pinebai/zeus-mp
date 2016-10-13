c=======================================================================
c/////////////////////////      RADIALPLOT      \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine radialplot 
c
c
c
c  FORMATTED WRITE OF theta-phi average (not integrated!) of SELECTED 
c  VARIABLES
c
c  PURPOSE: not only near the equator but full theta: from 0 to pi 
c  unlike fig. 6 of Stone & Pringle
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cosmic
      use brag
c
      implicit NONE
c
      integer, parameter :: nscal = 22 
      integer :: i, j, k, l
      real(rl) :: dvol, vol
      real(rl) :: scal(nscal)
      real(rl) :: dvr_av(in), pbyd_av(in), pcrbyd_av(in)
      real(rl) :: ddvr, dpbyd, dpcrbyd
      real(rl) :: b2cenx1, b3cenx1, bctr, bh1, bh2, bh3 
      real(rl) :: buf_send(3*in), buf_recv(3*in)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

      open(unit=41,file=oneDfile,status='unknown') 

      write(41,*) '#time = ', time, ' nhy =', nhy
      write(41,*)

      do i = is, ie+1
        dvr_av(i) = 0.0
        pbyd_av(i) = 0.0
        pcrbyd_av(i) = 0.0
        vol = 0.0        
      do k = ks, ke
      do j = js, je
        dvol = dvl2a(j)*dvl3a(k)
        vol = vol + dvol
c
c face centered
c
        dvr_av(i) = dvr_av(i) + dvol*0.5*(d(i,j,k)+d(i-1,j,k))*v1(i,j,k)
        pbyd_av(i) = pbyd_av(i) + dvol*(e(i,j,k)+e(i-1,j,k))*gamm1
     &  /(d(i,j,k)+d(i-1,j,k))
        if (xcosmic)   
     &   pcrbyd_av(i) = pcrbyd_av(i) + dvol*(ecr(i,j,k)+ecr(i-1,j,k))
     &  *gamm1_cr/(d(i,j,k)+d(i-1,j,k))
      enddo 
      enddo
      enddo
      dvr_av = dvr_av/vol; pbyd_av = pbyd_av/vol 
      pcrbyd_av = pcrbyd_av/vol

      buf_send(1:ie-is+2) = dvr_av(is:ie+1)
      buf_send(ie-is+3:2*ie-2*is+4) = pbyd_av(is:ie+1)
      buf_send(2*ie-2*is+5:3*ie-3*is+6) = pcrbyd_av(is:ie+1)    

#ifdef MPI_USED
      call MPI_ALLREDUCE( buf_send(1), buf_recv(1), 3*(ie-is+2)
     &, MPI_FLOAT, MPI_SUM, comm2d, ierr)
#endif

      dvr_av(is:ie+1) = buf_recv(1:ie-is+2)/(ntiles(2)*ntiles(3))
      pbyd_av(is:ie+1) = buf_recv(ie-is+3:2*ie-2*is+4)
     &/(ntiles(2)*ntiles(3))
      pcrbyd_av(is:ie+1) = buf_recv(2*ie-2*is+5:3*ie-3*is+6)
     &/(ntiles(2)*ntiles(3))

      do i = is, ie+1
        scal = 0.0d0
        vol = 0.0 
        do k=ks,ke
        do j=js,je
c
c quantities needed for flux calculation
c
              ddvr = 0.5*(d(i,j,k)+d(i-1,j,k))*v1(i,j,k) - dvr_av(i)
              dpbyd = (e(i,j,k)+e(i-1,j,k))*gamm1/(d(i,j,k)+d(i-1,j,k))
     &        - pbyd_av(i)
              if (xcosmic) dpcrbyd = (ecr(i,j,k)+ecr(i-1,j,k))*gamm1_cr
     &        /(d(i,j,k)+d(i-1,j,k)) - pcrbyd_av(i)
c 
              dvol = dvl2a(j)*dvl3a(k) 
c energy fluxes
c KE advection
              scal(11) = scal(11) + dvol*0.5*dvr_av(i)*(v1(i,j,k)**2 
     &  + 0.0625*(v2(i,j,k)+v2(i-1,j,k)+v2(i,j+1,k)+v2(i-1,j+1,k))**2 
     &  + 0.0625*(v3(i,j,k)+v3(i-1,j,k)+v3(i,j,k+1)+v3(i-1,j,k+1))**2)
c turbulent KE
              scal(12) = scal(12) + dvol*0.5*ddvr*(v1(i,j,k)**2
     &  + 0.0625*(v2(i,j,k)+v2(i-1,j,k)+v2(i,j+1,k)+v2(i-1,j+1,k))**2
     &  + 0.0625*(v3(i,j,k)+v3(i-1,j,k)+v3(i,j,k+1)+v3(i-1,j,k+1))**2)
c advective thermal enthalpy
              scal(13) = scal(13) + dvol*gamma/gamm1*pbyd_av(i)
     &        *dvr_av(i) 
c turbulent thermal enthalpy
              scal(14) = scal(14) + dvol*gamma/gamm1*dpbyd*ddvr
c advective CR enthalpy
              scal(15) = scal(15) + dvol*gamma_cr/gamm1_cr*pcrbyd_av(i)
     &        *dvr_av(i)
c turbulent CR enthalpy
              scal(16) = scal(16) + dvol*gamma_cr/gamm1_cr*dpcrbyd*ddvr
              
        if (xbrag) then  
          b2cenx1 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &    +b2(i-1,j,k)+b2(i-1,j+1,k))
          b3cenx1 = 0.25*(b3(i,j,k)+b3(i,j,k+1)
     &    +b3(i-1,j,k)+b3(i-1,j,k+1))
          bctr = sqrt(b1(i,j,k)*b1(i,j,k)
     &    +b2cenx1*b2cenx1+b3cenx1*b3cenx1)
          bh1 = b1(i,j,k)/(bctr+tiny)
          bh2 = b2cenx1/(bctr+tiny)
          bh3 = b3cenx1/(bctr+tiny)

c thermal conductive flux

              scal(17) = scal(17) - dvol*0.5*(chi(i,j,k)+chi(i-1,j,k))
     & *bh1*gamm1*( bh1*(e(i,j,k)/d(i,j,k)-e(i-1,j,k)/d(i-1,j,k))
     & *dx1bi(i)
     & +bh2*0.25*(e(i,j+1,k)/d(i,j+1,k)+e(i-1,j+1,k)/d(i-1,j+1,k)
     & -e(i,j-1,k)/d(i,j-1,k)-e(i-1,j-1,k)/d(i-1,j-1,k))
     & *g2ai(i)*dx2ai(j)
     & +bh3*0.25*(e(i,j,k+1)/d(i,j,k+1)+e(i-1,j,k+1)/d(i-1,j,k+1)
     & -e(i,j,k-1)/d(i,j,k-1)-e(i-1,j,k-1)/d(i-1,j,k-1))
     & *g31ai(i)*g32bi(j)*dx3ai(k) )

        endif

        if (xcosmic) then
 
         Dcr(i,j,k) = alpha_cr*x1a(i)*bctr
     &/sqrt(0.5*(d(i,j,k)+d(i-1,j,k)))


              scal(18) = scal(18) - dvol*Dcr(i,j,k)*bh1
     &   *( bh1*(ecr(i,j,k)-ecr(i-1,j,k))*dx1bi(i) 
     &   + bh2*0.25*(ecr(i,j+1,k)+ecr(i-1,j+1,k)
     &   -ecr(i,j-1,k)-ecr(i-1,j-1,k))*g2ai(i)*dx2ai(j)
     &   + bh3*0.25*(ecr(i,j,k+1)+ecr(i-1,j,k+1)
     &   -ecr(i,j,k-1)-ecr(i-1,j,k-1))*g31ai(i)*g32bi(j)*dx3ai(k) )
        endif
c |bhat.rhat|
        if (xbrag) then
              scal(19) =  scal(19) + dvol*abs(bh1)

c Spitzer value of heat flux
              scal(20) = scal(20) - dvol*0.5*(chi(i,j,k)+chi(i-1,j,k))
     &        *gamm1*(e(i,j,k)/d(i,j,k)-e(i-1,j,k)/d(i-1,j,k))*dx1bi(i)
        endif
c abundances
        if (nspec.gt.1) then
              scal(21) = scal(21) + dvol*abun(i,j,k,1)
              scal(22) = scal(22) + dvol*abun(i,j,k,2) 
        endif
c density
              scal(2) = scal(2) + dvol*d(i,j,k)
c internal energy density
              scal(3) = scal(3) + dvol*e(i,j,k)
              if (xmhd) then
c magnetic energy density
                scal(4) = scal(4) + dvol*0.5*( 0.5*(b1(i,j,k)**2
     &              +b1(i+1,j,k)**2+b2(i,j,k)**2+b2(i,j+1,k)**2)
     &              +b3(i,j,k)**2 )
c Maxwell stress
                scal(7) = scal(7) - dvol*0.5*(b1(i,j,k)+b1(i+1,j,k))
     &                *b3(i,j,k)
              endif 
c v3=r*sin(theta)*Omega
              scal(5) = scal(5) + dvol*v3(i,j,k)
c isothermal plasma sound speed
              scal(6) = scal(6) + dvol
     &                *sqrt((gamma-1.0)*e(i,j,k)/d(i,j,k)+tiny)
c radial velocity
              scal(8) = scal(8) + dvol*0.5*(v1(i,j,k)+v1(i+1,j,k))
c total Reynolds stress
              scal(9) = scal(9) + dvol*0.5*d(i,j,k)*(v1(i,j,k)
     &              +v1(i+1,j,k))*v3(i,j,k)
c cosmic ray energy density
              if (xcosmic) scal(10) = scal(10) + dvol*ecr(i,j,k)
              vol = vol + dvol
        enddo
        enddo
        scal = scal/(vol+tiny)
        scal(1) = x1b(i)
        write(41,2001) (scal(l), l = 1, nscal) 
      enddo
      close(41)
2001  format(57e20.7E4)
c
      return
      end
