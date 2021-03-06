c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   QUPDATE                 \\\\\\\\\\
c
c                            Developed by Prateek Sharma
c
c  called by driver qupdate_iso_d 
c 
c
c=======================================================================
c
      subroutine qupdate_iso(ibeg,iend,jbeg,jend,kbeg,kend,eold,enew
     &,dtcond)
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
      real(rl) :: bctr, q1(in,jn,kn), q2(in,jn,kn), q3(in,jn,kn)
      real(rl) :: divq, dTx1, dTx2, dTx3  
      real(rl) :: temp(in,jn,kn)
      real(rl) :: dtcond
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

          temp(i,j,k) = gamm1*eold(i,j,k)/d(i,j,k)
          chi(i,j,k) = kappa(i,j,k)*d(i,j,k)

        enddo
        enddo

      endif

      if (ldimen.eq.1) then

        k = ks
        j = js
        do i = ibeg-1, iend+1

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
         dTx1 = (temp(i,j,k)-temp(i-1,j,k))*dx1bi(i)
         q1(i,j,k) = -min(chi(i,j,k),chi(i-1,j,k))*dTx1
c
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
         dTx2 = (temp(i,j,k)-temp(i,j-1,k))*g2bi(i)*dx2bi(j)
         q2(i,j,k) = -min(chi(i,j,k),chi(i,j-1,k))*dTx2 
       enddo
       enddo
       enddo
c
c calculate q3
c
       do k = kbeg, kend+1
       do j = jbeg, jend
       do i = ibeg, iend
        dTx3 = (temp(i,j,k)-temp(i,j,k-1))*g31bi(i)*g32bi(j)*dx3bi(k)
        q3(i,j,k) = -min(chi(i,j,k),chi(i,j,k-1))*dTx3
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
