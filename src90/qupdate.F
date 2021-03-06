c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   QUPDATE                 \\\\\\\\\\
c
c                            Developed by Prateek Sharma
c
c  called by driver qupdate_d 
c 
c
c=======================================================================
c
      subroutine qupdate(ibeg,iend,jbeg,jend,kbeg,kend,eold,enew,dtcond)
c
c......................................................................
      use real_prec
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
      use config
      use brag
c
      implicit NONE
c
      integer  :: i, j, k
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend 
      real(rl) :: bh1(in,jn,kn), bh2(in,jn, kn), bh3(in,jn,kn)   
      real(rl) :: bh2cenx1(in,jn,kn), bh3cenx1(in,jn,kn)
      real(rl) :: bh1cenx2(in,jn,kn), bh3cenx2(in,jn,kn)
      real(rl) :: bh1cenx3(in,jn,kn), bh2cenx3(in,jn,kn) 
      real(rl) :: b1cenx2, b1cenx3, b2cenx1, b2cenx3, b3cenx1, b3cenx2
      real(rl) :: bctr, q1(in,jn,kn), q2(in,jn,kn), q3(in,jn,kn)
      real(rl) :: divq  
      real(rl) :: temp(in,jn,kn)
      real(rl) :: dtcond
      real(rl) :: dTx1, dTx2, dTx3
      real(rl) :: lim 
      real(rl) :: eold(in,jn,kn), enew(in,jn,kn)
c
c----------------------------------------------------------------------
c compute magnetic field unit vectors
c
c----------------------------------------------------------------------
c

      if (ldimen.eq.3) then

        do k = kbeg-1, kend+1
        do j = jbeg-1, jend+1
        do i = ibeg-1, iend+1
          b2cenx1 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &    +b2(i-1,j,k)+b2(i-1,j+1,k))
          b3cenx1 = 0.25*(b3(i,j,k)+b3(i,j,k+1)
     &    +b3(i-1,j,k)+b3(i-1,j,k+1))
          bctr = sqrt(b1(i,j,k)*b1(i,j,k)
     &    +b2cenx1*b2cenx1+b3cenx1*b3cenx1)
          bh1(i,j,k) = b1(i,j,k)/(bctr+tiny)
          bh2cenx1(i,j,k) = b2cenx1/(bctr+tiny)
          bh3cenx1(i,j,k) = b3cenx1/(bctr+tiny) 
          b1cenx2 = 0.25*(b1(i,j,k)+b1(i+1,j,k)
     &    +b1(i,j-1,k)+b1(i+1,j-1,k))
          b3cenx2 = 0.25*(b3(i,j,k)+b3(i,j,k+1)
     &    +b3(i,j-1,k)+b3(i,j-1,k+1))
          bctr = sqrt(b1cenx2*b1cenx2+b2(i,j,k)*b2(i,j,k)
     &    +b3cenx2*b3cenx2)
          bh1cenx2(i,j,k) = b1cenx2/(bctr+tiny)
          bh2(i,j,k) = b2(i,j,k)/(bctr+tiny)
          bh3cenx2(i,j,k) = b3cenx2/(bctr+tiny)
          b1cenx3 = 0.25*(b1(i,j,k)+b1(i+1,j,k)
     &    +b1(i,j,k-1)+b1(i+1,j,k-1))
          b2cenx3 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &    +b2(i,j,k-1)+b2(i,j+1,k-1))
          bctr = sqrt(b1cenx3*b1cenx3+b2cenx3*b2cenx3
     &    +b3(i,j,k)*b3(i,j,k))
          bh1cenx3(i,j,k) = b1cenx3/(bctr+tiny)
          bh2cenx3(i,j,k) = b2cenx3/(bctr+tiny)
          bh3(i,j,k) = b3(i,j,k)/(bctr+tiny)

          temp(i,j,k) = gamm1*eold(i,j,k)/d(i,j,k)
          chi(i,j,k) = kappa(i,j,k)*d(i,j,k)

        enddo
        enddo
        enddo
     
      endif


      if (ldimen.eq.2) then

        k = ks 
        do j = jbeg-1, jend+1
        do i = ibeg-1, iend+1
          b2cenx1 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &    +b2(i-1,j,k)+b2(i-1,j+1,k))
          b3cenx1 = 0.5*(b3(i,j,k)+b3(i-1,j,k))
          bctr = sqrt(b1(i,j,k)*b1(i,j,k)
     &    +b2cenx1*b2cenx1+b3cenx1*b3cenx1)
          bh1(i,j,k) = b1(i,j,k)/(bctr+tiny)
          bh2cenx1(i,j,k) = b2cenx1/(bctr+tiny)
          bh3cenx1(i,j,k) = b3cenx1/(bctr+tiny)
          b1cenx2 = 0.25*(b1(i,j,k)+b1(i+1,j,k)
     &    +b1(i,j-1,k)+b1(i+1,j-1,k))
          b3cenx2 = 0.5*(b3(i,j,k)+b3(i,j-1,k))
          bctr = sqrt(b1cenx2*b1cenx2+b2(i,j,k)*b2(i,j,k)
     &    +b3cenx2*b3cenx2)
          bh1cenx2(i,j,k) = b1cenx2/(bctr+tiny)
          bh2(i,j,k) = b2(i,j,k)/(bctr+tiny)
          bh3cenx2(i,j,k) = b3cenx2/(bctr+tiny)
          b1cenx3 = 0.5*(b1(i,j,k)+b1(i+1,j,k))
          b2cenx3 = 0.5*(b2(i,j,k)+b2(i,j+1,k))
          bctr = sqrt(b1cenx3*b1cenx3+b2cenx3*b2cenx3
     &    +b3(i,j,k)*b3(i,j,k))
          bh1cenx3(i,j,k) = b1cenx3/(bctr+tiny)
          bh2cenx3(i,j,k) = b2cenx3/(bctr+tiny)
          bh3(i,j,k) = b3(i,j,k)/(bctr+tiny)

          temp(i,j,k) = gamm1*eold(i,j,k)/d(i,j,k)
          chi(i,j,k) = kappa(i,j,k)*d(i,j,k)

        enddo
        enddo

      endif

      if (ldimen.eq.1) then

        k = ks
        j = js
        do i = ibeg-1, iend+1
          b2cenx1 = 0.5*(b2(i,j,k)+b2(i-1,j,k))
          b3cenx1 = 0.5*(b3(i,j,k)+b3(i-1,j,k))
          bctr = sqrt(b1(i,j,k)*b1(i,j,k)+b2cenx1*b2cenx1
     &    +b3cenx1*b3cenx1)
          bh1(i,j,k) = b1(i,j,k)/(bctr+tiny)
          bh2cenx1(i,j,k) = b2cenx1/(bctr+tiny)
          bh3cenx1(i,j,k) = b3cenx1/(bctr+tiny)
          b1cenx2 = 0.5*(b1(i,j,k)+b1(i+1,j,k))
          b3cenx2 = b3(i,j,k)
          bctr = sqrt(b1cenx2*b1cenx2+b2(i,j,k)*b2(i,j,k)
     &    +b3cenx2*b3cenx2)
          bh1cenx2(i,j,k) = b1cenx2/(bctr+tiny)
          bh2(i,j,k) = b2(i,j,k)/(bctr+tiny)
          bh3cenx2(i,j,k) = b3cenx2/(bctr+tiny)
          b1cenx3 = 0.5*(b1(i,j,k)+b1(i+1,j,k))
          b2cenx3 = b2(i,j,k)
          bctr = sqrt(b1cenx3*b1cenx3+b2cenx3*b2cenx3
     &    +b3(i,j,k)*b3(i,j,k))
          bh1cenx3(i,j,k) = b1cenx3/(bctr+tiny)
          bh2cenx3(i,j,k) = b2cenx3/(bctr+tiny)
          bh3(i,j,k) = b3(i,j,k)/(bctr+tiny)

          temp(i,j,k) = gamm1*eold(i,j,k)/d(i,j,k)
          chi(i,j,k) = kappa(i,j,k)*d(i,j,k)

        enddo

      endif

c
c calculate q1
c
       do k = kbeg, kend
       do j = jbeg, jend
       do i = ibeg, iend+1
c
c done slightly differently [less computationally intensive] than
c in the axisymmetric code
c
         dTx2 = lim( lim(temp(i,j+1,k)-temp(i,j,k)
     &             ,temp(i-1,j+1,k)-temp(i-1,j,k))
     &             , lim(temp(i,j,k)-temp(i,j-1,k)
     &             ,temp(i-1,j,k)-temp(i-1,j-1,k)) )
         dTx3 = lim( lim(temp(i,j,k+1)-temp(i,j,k)
     &             ,temp(i-1,j,k+1)-temp(i-1,j,k))
     &             , lim(temp(i,j,k)-temp(i,j,k-1)
     &             ,temp(i-1,j,k)-temp(i-1,j,k-1)) )
         dTx1 = (temp(i,j,k)-temp(i-1,j,k))*dx1bi(i)
         dTx2 = dTx2*g2ai(i)*dx2ai(j)
         dTx3 = dTx3*g31ai(i)*g32bi(j)*dx3ai(k)
         q1(i,j,k) = -min(chi(i,j,k),chi(i-1,j,k))*bh1(i,j,k)
     &   *( bh1(i,j,k)*dTx1 + bh2cenx1(i,j,k)*dTx2 
     &   + bh3cenx1(i,j,k)*dTx3 )
       enddo
       enddo
       enddo 
c
c calculate q2
c       
       do k = kbeg, kend
c
c keeping a spherical polar coordinate system in mind; avoid singularities at poles
c
       do j = jbeg, jend+1
       do i = ibeg, iend
         dTx1 = lim( lim(temp(i+1,j,k)-temp(i,j,k)
     &          ,temp(i+1,j-1,k)-temp(i,j-1,k))
     &          , lim(temp(i,j,k)-temp(i-1,j,k)
     &          ,temp(i,j-1,k)-temp(i-1,j-1,k)) )
         dTx3 = lim( lim(temp(i,j,k+1)-temp(i,j,k)
     &          ,temp(i,j-1,k+1)-temp(i,j-1,k))
     &          , lim(temp(i,j,k)-temp(i,j,k-1)
     &          ,temp(i,j-1,k)-temp(i,j-1,k-1)) )
         dTx1 = dTx1*dx1ai(i)
         dTx2 = (temp(i,j,k)-temp(i,j-1,k))*g2bi(i)*dx2bi(j)
         dTx3 = dTx3*g31bi(i)*g32ai(j)*dx3ai(k)
         if (x2a(j).le.1d-10.and.x2a(j).ge.-1d-10) dTx3 = 0.0d0 
         if (x2a(j).le.pi+1d-10.and.x2a(j).ge.pi-1d-10) dTx3 = 0.0d0
         q2(i,j,k) = -min(chi(i,j,k),chi(i,j-1,k))*bh2(i,j,k)
     &   *( bh1cenx2(i,j,k)*dTx1 + bh2(i,j,k)*dTx2 
     &   + bh3cenx2(i,j,k)*dTx3 ) 
       enddo
       enddo
       enddo
c
c calculate q3
c
       do k = kbeg, kend+1
       do j = jbeg, jend
       do i = ibeg, iend
         dTx1 = lim( lim(temp(i+1,j,k)-temp(i,j,k)
     &             ,temp(i+1,j,k-1)-temp(i,j,k-1))
     &             , lim(temp(i,j,k)-temp(i-1,j,k)
     &             ,temp(i,j,k-1)-temp(i-1,j,k-1)) )
         dTx2 = lim( lim(temp(i,j+1,k)-temp(i,j,k)
     &          ,temp(i,j+1,k-1)-temp(i,j,k-1))
     &          , lim(temp(i,j,k)-temp(i,j-1,k)
     &          ,temp(i,j,k-1)-temp(i,j-1,k-1)) )
        dTx1 = dTx1*dx1ai(i) 
        dTx2 = dTx2*g2bi(i)*dx2ai(j)
        dTx3 = (temp(i,j,k)-temp(i,j,k-1))*g31bi(i)*g32bi(j)*dx3bi(k)
        q3(i,j,k) = -min(chi(i,j,k),chi(i,j,k-1))*bh3(i,j,k)
     &  *( bh1cenx3(i,j,k)*dTx1 + bh2cenx3(i,j,k)*dTx2
     &  + bh3(i,j,k)*dTx3 )
       enddo
       enddo
       enddo
c
c update internal energy
c
       do k = kbeg, kend
       do j = jbeg, jend
       do i = ibeg, iend
         divq = (g2a(i+1)*g31a(i+1)*q1(i+1,j,k)
     &   -g2a(i)*g31a(i)*q1(i,j,k))*dvl1ai(i)
     &   + (g32a(j+1)*q2(i,j+1,k)-g32a(j)*q2(i,j,k))*dvl2ai(j)*g2bi(i)
     &   + (q3(i,j,k+1)-q3(i,j,k))*dvl3ai(k)*g31bi(i)*g32bi(j)
         enew(i,j,k) = eold(i,j,k) - dtcond*divq
       enddo
       enddo
       enddo
c
      return
      end

c
c MC limiter is used as lim {it worked best in our anisotropic conduction 
c paper}
c 

      function lim(var1, var2)
      use real_prec 
      implicit NONE
      real(rl) :: var1, var2, temp, lim, sgn, absv1, absv2

      if (var1*var2 .le. 0.0) then
        lim = 0.0d0
        return       
      else
        sgn = sign(1.0d0,var1)
        absv1 = abs(var1)
        absv2 = abs(var2)
        temp = sgn*2.0*min(absv1,absv2)
      endif

      lim = sgn*min(abs(temp), 0.5*(absv1+absv2))

      end function lim
