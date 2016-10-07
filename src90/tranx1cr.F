c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                TRANX1CR                   \\\\\\\\\c
c                            Developed by
c                           Prateek Sharma
c
c=======================================================================
c
      subroutine tranx1cr (ibeg,iend,jbeg,jend,kbeg,kend,dtstr)
c
      use real_prec
      use param
      use config
      use field
      use grid
      use mpiyes
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: dqm, dqp
      real(rl) :: atwid(ijkn),
     .            ecrst_twid(ijkn), decrst(ijkn), xi_st, q1_st, vst
      real(rl) :: xarg
      real(rl) :: sgn, dtstr, bdotpcrx1, bdotpcrx2, bdotpcrx3
c-----------------------------------------------------------------------
c
      do 10 i=ibeg-1,iend+1
       if(xvgrid) then
        atwid (i)  =       g2ah(i) * g31ah(i)
       else
        atwid (i)  =       g2a (i) * g31a (i)
       endif
10    continue
c
      do 100 k=kbeg,kend
       do 90 j=jbeg,jend
c
        if (ibeg .eq. is) then  !  Need d(is-3) from neighbor.
         i         = is - 2
         decrst(i) = zro
        endif ! ibeg
        do 30 i=max(ibeg-2,is-1),iend+1
          dqm       = (ecr(i  ,j,k) - ecr(i-1,j,k)) * dx1bi(i  )
          dqp       = (ecr(i+1,j,k) - ecr(i  ,j,k)) * dx1bi(i+1)
          decrst(i) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
30      continue
c
        do 40 i=ibeg-1,iend+1
          bdotpcrx1 = b1(i,j,k)*(ecr(i,j,k)-ecr(i-1,j,k))*dx1bi(i)
          bdotpcrx2 = 0.25*b2(i,j,k)*(ecr(i,j,k)-ecr(i,j-1,k))*dx2bi(j)
     &    *g2bi(i)
     &    + 0.25*b2(i,j+1,k)*(ecr(i,j+1,k)-ecr(i,j,k))*dx2bi(j+1)
     &    *g2bi(i)
     &    + 0.25*b2(i-1,j,k)*(ecr(i-1,j,k)-ecr(i-1,j-1,k))*dx2bi(j)
     &    *g2bi(i-1)
     &    + 0.25*b2(i-1,j+1,k)*(ecr(i-1,j+1,k)-ecr(i-1,j,k))*dx2bi(j+1)
     &    *g2bi(i-1)
          bdotpcrx3 = 0.25*b3(i,j,k)*(ecr(i,j,k)-ecr(i,j,k-1))*dx3bi(k)
     &    *g31bi(i)*g32bi(j)
     &    + 0.25*b3(i,j,k+1)*(ecr(i,j,k+1)-ecr(i,j,k))*dx3bi(k+1)
     &    *g31bi(i)*g32bi(j)
     &    + 0.25*b3(i-1,j,k)*(ecr(i-1,j,k)-ecr(i-1,j,k-1))*dx3bi(k)
     &    *g31bi(i-1)*g32bi(j)
     &    + 0.25*b3(i-1,j,k+1)*(ecr(i-1,j,k+1)-ecr(i-1,j,k))*dx3bi(k+1)
     &    *g31bi(i-1)*g32bi(j) 
c          xarg = (-bdotpcrx1-bdotpcrx2-bdotpcrx3)/(0.1*1.0*8.0)
c          vst = v1(i,j,k) - vg1(i)
c     &    + tanh(xarg)*b1(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i-1,j,k))) 
          vst = v1(i,j,k) - vg1(i)
     &    + sgn( -bdotpcrx1-bdotpcrx2-bdotpcrx3 )
     &    *b1(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i-1,j,k)))
          xi_st     =  dtstr*vst 
          q1_st     = haf*sgn( xi_st )
          ecrst_twid (i) = ( 0.5 + q1_st ) * ( ecr(i-1,j,k)
     1              + ( dx1a(i-1) - xi_st ) * decrst (i-1) )
     2              + ( 0.5 - q1_st ) * ( ecr(i  ,j,k)
     3              - ( dx1a(i  ) + xi_st ) * decrst (i  ) )
          ecrst_twid(i) = xi_st*atwid(i)*ecrst_twid(i)
40      continue
c
        do 80 i=ibeg,iend
          if (xvgrid) then
             ecr(i,j,k) = ( ecr  (i,j,k) * dvl1a(i)
     1               - ecrst_twid(i+1) + ecrst_twid (i) )* dvl1ani(i)
          else
             ecr(i,j,k) = ( ecr  (i,j,k) * dvl1a(i)
     1               - ecrst_twid(i+1) + ecrst_twid (i) )* dvl1a i(i)
          endif
80        continue
90      continue
100    continue
c       
      return
      end
c
      function sgn(var)
      use param
      implicit NONE
      real*8 :: sgn, var
       if (var.gt.10.*tiny) sgn=1.0
       if (10.*tiny.ge.var.and.var.le.10.*tiny) sgn=0.0
       if (var.lt.-tiny*10.) sgn=-1.0
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                TRANX1CR                   \\\\\\\\\c
c=======================================================================
