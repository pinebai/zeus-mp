c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 SHOCK_TUBE                \\\\\\\\\\
c
c=======================================================================
c
      subroutine kh
c
c
c written by Prateek Sharma
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use cons
      use metallicity
      use gravmod
      use root
      use domain
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit none
c
      integer  :: i, j, k
      real(rl) :: amp
      namelist  /pgen/ amp 
c
c-----------------------------------------------------------------------
c
      do i = 1, jn

c#ifdef 0
c        if (x2b(i).lt.0.25) then
c         	d(:,i,:) = 1.0
c         	p(:,i,:) = 2.5
c         	pcr(i,:,:) = 1.3e5
c         	abun(i,:,:,1) = 1.0
c        else
c        	d(:,i,:) = 2.0
c        	p(:,i,:) = 2.5
c       	pcr(:,i,:) = 2.4e2
c        	abun(i,:,:,1) = 1.e-15
c        endif  
c        v1(i,:,:) = 0.0
c#endif
	if (x2b(i).lt.0.25 .or. x2b(i).gt. 0.75) then
		d(:,i,:) = 1.0
		p(:,i,:) = 2.5
c		pcr(i,:,:) = 1.3e5
c		abun(i,:,:,1) = 1.0
		v1(:,i,:) = 0.5
	else
		d(:,i,:) = 2.0
		p(:,i,:) = 2.5
		v1(:,i,:) = -0.5
c		pcr(:,i,:) = 2.4e2
c		abun(i,:,:,1) = 1.0
	endif


c        if (x1b(i).lt.0.5) then
c         d(i,:,:) = 1.0
c          p(i,:,:) = 1.0/2.
c         pcr(i,:,:) = 1.0/2.
c          abun(i,:,:,1) = 1.0
c        else
c          d(i,:,:) = 0.125
c          p(i,:,:) = 0.1/2.
c         pcr(i,:,:) = 0.1/2.
c          abun(i,:,:,1) = 1.e-15
c        endif
c        v1(:,i,:) =0.0

      enddo
	do i=1,in
		v2(i,:,:) = 0.0001 * sin(4*pi*x1b(i))
	enddo

      e = p/gamm1
c      ecr = pcr/gamm1_cr
c
      return
      end
