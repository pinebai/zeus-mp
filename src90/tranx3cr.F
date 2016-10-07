crc=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                TRANX3CR                   \\\\\\\\\c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine tranx3cr (ibeg,iend,jbeg,jend,kbeg,kend,dtstr)
c
      use real_prec
      use config
      use param
      use field
      use grid
      use mpiyes
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: dqm, dqp, xi_st, q1_st
c
      real(rl) :: atwid2(ijkn),  atwid(ijkn)
      real(rl) :: ecrst_twid (ijkn,1), decrst(ijkn,1)
      real(rl) :: sgn, dtstr, vst, bdotpcrx1, bdotpcrx2, bdotpcrx3
c-----------------------------------------------------------------------
c
      do 10 i=ibeg-1,iend+1
       atwid (i)           = g2b(i) * dx1a(i) * dvl1ai(i)
10    continue
      do 20 j=jbeg-1,jend+1
       atwid2(j)           = dx2a(j) * dvl2ai(j)
20    continue
c
      do 2100 j=jbeg,jend
         do 1090 i=ibeg,iend
c
          if (kbeg .eq. ks) then  !  Need d(ks-3) from neighbor.
           k        = ks - 2
           decrst(k,1) = zro
          endif
          do 1030 k=max(kbeg-2,ks-1),kend+1
            dqm      = (ecr(i  ,j,k  ) - ecr(i  ,j,k-1))* dx3bi(k  )
            dqp      = (ecr(i  ,j,k+1) - ecr(i  ,j,k  ))* dx3bi(k+1)
            decrst(k,1)= max ( dqm * dqp, zro )
     1               * sign ( one, dqm + dqp )
     2               / max ( abs ( dqm + dqp ), tiny )
1030      continue
c
          do 1040 k=kbeg-1,kend+1
           bdotpcrx1 = 0.25*b1(i,j,k)*(ecr(i,j,k)-ecr(i-1,j,k))*dx1bi(i)
     &     + 0.25*b1(i+1,j,k)*(ecr(i+1,j,k)-ecr(i,j,k))*dx1bi(i+1)
     &     + 0.25*b1(i,j,k-1)*(ecr(i,j,k-1)-ecr(i-1,j,k-1))*dx1bi(i)
     &     + 0.25*b1(i+1,j,k-1)*(ecr(i+1,j,k-1)-ecr(i,j,k-1))*dx1bi(i+1)           
           bdotpcrx2 = 0.25*b2(i,j,k)*(ecr(i,j,k)-ecr(i,j-1,k))*dx2bi(j)
     &     *g2bi(i)
     &     + 0.25*b2(i,j+1,k)*(ecr(i,j+1,k)-ecr(i,j,k))*dx2bi(j+1)
     &     *g2bi(i)
     &     + 0.25*b2(i,j,k-1)*(ecr(i,j,k-1)-ecr(i,j-1,k-1))*dx2bi(j)
     &     *g2bi(i)
     &     + 0.25*b2(i,j+1,k-1)*(ecr(i,j+1,k-1)-ecr(i,j,k-1))*dx2bi(j+1)
     &     *g2bi(i)
           bdotpcrx3 = b3(i,j,k)*(ecr(i,j,k)-ecr(i,j,k-1))*dx3bi(k)
     &     *g31bi(i)*g32bi(j)
           
           vst = v3(i,j,k) - vg3(k) 
     &    + sgn( -b3(i,j,k)*(ecr(i,j,k)-ecr(i,j,k-1)) )
     &    *b3(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i,j,k-1)))
           xi_st = dtstr*vst
           q1_st = sgn ( haf, xi_st )
            ecrst_twid (k,1) = ( 0.5 + q1_st ) * ( ecr(i  ,j,k-1)
     1                  + ( dx3a(k-1) - xi_st ) * decrst (k-1,1) )
     2                  + ( 0.5 - q1_st ) * ( ecr(i  ,j,k  )
     3                  - ( dx3a(k  ) + xi_st ) * decrst (k  ,1) )
c
            ecrst_twid (k,1) = xi_st*atwid(i)*atwid2(j)*ecrst_twid (k,1)
1040      continue
c
           do 1080 k=kbeg,kend
               if(xvgrid) then
                ecr(i,j,k) = ( ecr  (i  ,j,k) * dvl3a(k)
     .          - ecrst_twid (k+1,1) + ecrst_twid (k,1) ) * dvl3ani(k)
               else
                ecr(i,j,k) = ( ecr  (i  ,j,k) * dvl3a(k)
     .          - ecrst_twid (k+1,1) + ecrst_twid (k,1) ) * dvl3a i(k)
               endif
1080      continue
1090    continue
2100  continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                T R A N X 3                \\\\\\\\\c
c=======================================================================
c
c
