c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                TRANX2CR_D                 \\\\\\\\\c
c                            Developed by
c                           Prateek Sharma
c
c=======================================================================
c
       subroutine tranx2cr_d
c
c driver for tranx2cr which advects cr at va (modulo grad(pcr))
c
c-----------------------------------------------------------------------
c
      use real_prec
      use param
      use root
      use field
      use grid
      use bndry
      use mpiyes
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, l, ncr
c
      real(rl) :: decrp, decrm, Lcr, dtstr, valf
c-----------------------------------------------------------------------
c
      dtstr = huge 
      do k = ks, ke
      do j = js, je
      do i = is, ie
        decrp = abs(ecr(i,j+1,k)-ecr(i,j,k))*dx2bi(j+1)*g2bi(i)
        decrm = abs(ecr(i,j,k)-ecr(i,j-1,k))*dx2bi(j)*g2bi(i)
        valf = max( abs(b2(i,j,k))/sqrt(0.5*(d(i,j,k)+d(i,j-1,k))),
     &  abs(b2(i,j+1,k))/sqrt(0.5*(d(i,j,k)+d(i,j+1,k))) )
c taking a upper bound Lcr to be 200kpc; must experiment with this
        Lcr = min(ecr(i,j,k)/max(decrm,decrp), 6.173e23)
c        Lcr = min(ecr(i,j,k)/max(decrm,decrp), 20.0)
        dtstr = min(dtstr, (g2b(i)*dx2a(j))**2/(valf*Lcr))
c        write(*,201) decrp, decrm, valf, Lcr, dtstr
      enddo
      enddo
      enddo
201   format(5e20.7)
c
#ifdef MPI_USED
      buf_in(1) = dtstr
      call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
      dtstr  =   buf_out(1)
#endif
c
      ncr = int(dt/dtstr)+1
      ncr = min(ncr, 1000)
      dtstr = dt/float(ncr)
c      write(22,*) time, ncr
      nreq = 0
      nsub = nsub + 1   
      call bvalecr(1,1,3,3,1,1,ecr)
        if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
      do l = 1, ncr
        nreq = 0
        nsub = nsub + 1
        call tranx2cr(is, ie, js, je, ks, ke, dtstr)
        bvstat(1:7,9) = 0 ! ecr boundary out of date
c        call bvalecrstr(0,0,3,3,0,0,ecr)
        call bvalecrstr(1,1,3,3,1,1,ecr)
        if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
      enddo
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                TRANX2CR_D                \\\\\\\\\c
c=======================================================================
c
c
