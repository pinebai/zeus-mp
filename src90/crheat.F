c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 CRHEAT                    \\\\\\\\\c
c                            Developed by Prateek Sharma
c
c=======================================================================
c
      subroutine crheat
c
c......................................................................
      use real_prec
      use param
      use root
      use grid
      use field
      use mpiyes
      use mpipar
      use config
      use cosmic 
      use cons
      use bndry
c
      implicit NONE
      integer  :: i, j, k
      real(rl) :: heating_rate, valfp, valfm, dpcrdrp, dpcrdrm
c density has been updated in qupdate_d, and crdiffuse_d!      
      do k = ks, ke
      do j = js, je
      do i = is, ie
        valfm = b1(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i-1,j,k)))
        valfp = b1(i+1,j,k)/sqrt(0.5*(d(i,j,k)+d(i+1,j,k)))
        dpcrdrm = gamm1_cr*(ecr(i,j,k)-ecr(i-1,j,k))*dx1bi(i)
        dpcrdrp = gamm1_cr*(ecr(i+1,j,k)-ecr(i,j,k))*dx1bi(i+1)
        heating_rate = 0.5*valfm*dpcrdrm + 0.5*valfp*dpcrdrp
        valfm = b2(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i,j-1,k)))
        valfp = b2(i,j+1,k)/sqrt(0.5*(d(i,j,k)+d(i,j+1,k)))
        dpcrdrm = gamm1_cr*(ecr(i,j,k)-ecr(i,j-1,k))*dx2bi(j)*g2bi(i)
        dpcrdrp = gamm1_cr*(ecr(i,j+1,k)-ecr(i,j,k))*dx2bi(j+1)*g2bi(i)
        heating_rate = heating_rate + 0.5*valfm*dpcrdrm 
     &  + 0.5*valfp*dpcrdrp
        valfm = b3(i,j,k)/sqrt(0.5*(d(i,j,k)+d(i,j,k-1)))
        valfp = b3(i,j,k+1)/sqrt(0.5*(d(i,j,k)+d(i,j,k+1)))
        dpcrdrm = gamm1_cr*(ecr(i,j,k)-ecr(i,j,k-1))*dx3bi(k)*g31bi(i)
     &  *g32bi(j)
        dpcrdrp = gamm1_cr*(ecr(i,j,k+1)-ecr(i,j,k))*dx3bi(k+1)*g31bi(i)
     &  *g32bi(j)
        heating_rate =  abs(heating_rate + 0.5*valfm*dpcrdrm 
     &  + 0.5*valfp*dpcrdrp)
        e(i,j,k) = e(i,j,k) + dt*heating_rate
        ecr(i,j,k) = ecr(i,j,k) - dt*heating_rate
      enddo
      enddo
      enddo
c
c making the boundary values out of date
c
      do i = 1,6
       bvstat(i,2) = 0      !  e
       bvstat(i,9) = 0      !  ecr
      enddo
c
      return
      end
