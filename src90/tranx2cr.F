c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                TRANX2CR                   \\\\\\\\\c
c                            Developed by
c                           Prateek Sharma
c
c=======================================================================
c
       subroutine tranx2cr (ibeg,iend,jbeg,jend,kbeg,kend,dtstr)
c
      use real_prec
      use config
      use param
      use field
      use grid
      use root
      use mpiyes
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: dqm,dqp, xi_st,q1_st, atwid(ijkn)
c
      real(rl) :: ecrst_twid(ijkn), vst, decrst(ijkn)
      real(rl) :: sgn, dtstr, bdotpcrx1, bdotpcrx2, bdotpcrx3 
c
c-----------------------------------------------------------------------
c
c Compute time-centered area factors.
c
      do 10 i=ibeg-1,iend+1
        atwid (i)         = g31b(i) * dx1a(i) * dvl1ai(i)
10    continue
c
      do 2100 k=kbeg,kend
       do 1090 i=ibeg,iend
c
        if (jbeg .eq. js) then  !  Need d(js-3) from neighbor.
         j         = js - 2
         decrst(j) = zro
        endif
        do 1030 j=max(jbeg-2,js-1),jend+1
          dqm       = (ecr(i,j  ,k) - ecr(i,j-1,k)) * dx2bi(j  )
          dqp       = (ecr(i,j+1,k) - ecr(i,j  ,k)) * dx2bi(j+1)
          decrst(j) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
c          write(*,*) nhy, decrst(j)
1030    continue
c
        do 1040 j=jbeg-1,jend+1
          bdotpcrx1 = 0.25*b1(i,j,k)*(ecr(i,j,k)-ecr(i-1,j,k))*dx1bi(i)
     &    + 0.25*b1(i+1,j,k)*(ecr(i+1,j,k)-ecr(i,j,k))*dx1bi(i+1)
     &    + 0.25*b1(i,j-1,k)*(ecr(i,j-1,k)-ecr(i-1,j-1,k))*dx1bi(i)
     &    + 0.25*b1(i+1,j-1,k)*(ecr(i+1,j-1,k)-ecr(i,j-1,k))*dx1bi(i+1) 
          bdotpcrx2 = b2(i,j,k)*(ecr(i,j,k)-ecr(i,j-1,k))*dx2bi(j)
     &    *g2bi(i)
          bdotpcrx3 = 0.25*b3(i,j,k)*(ecr(i,j,k)-ecr(i,j,k-1))*dx3bi(k)
     &    *g31bi(i)*g32bi(j)
     &    + 0.25*b3(i,j,k+1)*(ecr(i,j,k+1)-ecr(i,j,k))*dx3bi(k+1)
     &    *g31bi(i)*g32bi(j)
     &    + 0.25*b3(i,j-1,k)*(ecr(i,j-1,k)-ecr(i,j-1,k-1))*dx3bi(k)
     &    *g31bi(i)*g32bi(j-1)
     &    + 0.25*b3(i,j-1,k+1)*(ecr(i,j-1,k+1)-ecr(i,j-1,k))*dx3bi(k+1)
     &    *g31bi(i)*g32bi(j-1)
c          write(*,*) nhy, bdotpcrx1, bdotpcrx2, bdotpcrx3
          vst = v2(i,j,k) - vg2(j)
     &    + sgn( -bdotpcrx1-bdotpcrx2-bdotpcrx3 )
     &    *b2(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i,j-1,k)))
          xi_st     =  dtstr*vst
          q1_st     = haf*sgn( xi_st )
          ecrst_twid (j) = ( 0.5 + q1_st ) * ( ecr(i,j-1,k)
     1              + ( dx2a(j-1) - xi_st ) * decrst (j-1) )
     2              + ( 0.5 - q1_st ) * ( ecr(i  ,j,k)
     3              - ( dx2a(j  ) + xi_st ) * decrst (j  ) )
          if(xvgrid) then
            ecrst_twid(j) = xi_st*g32ah(j)*atwid(j)*ecrst_twid(j)
          else
            ecrst_twid(j) = xi_st*g32a(j)*atwid(j)*ecrst_twid(j)
          endif
c          write(*,*) nhy, ecrst_twid(j) 
1040    continue
c
         do 1080 j=jbeg,jend
c
           if(xvgrid) then
            ecr  (i,j,k) = ( ecr  (i,j,k) * dvl2a(j)
     .             - ecrst_twid (j+1) + ecrst_twid (j)  ) * dvl2ani(j)
           else
            ecr  (i,j,k) = ( ecr  (i,j,k) * dvl2a(j)
     .           - ecrst_twid (j+1)  + ecrst_twid (j)  ) * dvl2a i(j)
           endif
1080     continue
1090   continue
2100   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                T R A N X 2                \\\\\\\\\c
c=======================================================================
c
c
