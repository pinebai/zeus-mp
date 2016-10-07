c=======================================================================
c
      subroutine floor_d 
c
c  PURPOSE: applies density floor  
c
c-----------------------------------------------------------------------
      use field
      use grid
      use root
c
      implicit NONE
c
      integer :: i, j, k
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

      do k = ks-2, ke+3
      do j = js-2, je+3
      do i = is-2, ie+3

        if (d(i,j,k).lt.dfloor) then
c
c madd is the mass added in real zones (not including the ghost zones)
c
          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je 
     &        .and.k.ge.ks.and.k.le.ke)
     &      madd = madd + (dfloor-d(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k) 
c
c but the floor is applied everywhere, including the ghosts
c
          d(i,j,k) = dfloor
        endif

      enddo 
      enddo 
      enddo
c
      return
      end

c=======================================================================
c
      subroutine floor_e
c
c  PURPOSE: applies internal energy floor
c
c-----------------------------------------------------------------------
      use field
      use grid
      use root
c
      implicit NONE
c
      integer :: i, j, k
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

      do k = ks-2, ke+3
      do j = js-2, je+3
      do i = is-2, ie+3

        if (e(i,j,k).lt.efloor) then
c
c madd is the mass added in real zones (not including the ghost zones)
c
          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je
     &        .and.k.ge.ks.and.k.le.ke)
     &      eadd = eadd + (efloor-e(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k)
c
c but the floor is applied everywhere, including the ghosts
c
          e(i,j,k) = efloor
        endif

      enddo
      enddo
      enddo
c
      return
      end
