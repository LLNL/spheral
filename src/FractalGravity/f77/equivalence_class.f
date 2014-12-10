      subroutine equivalence_class(nf,n,lista,listb,m,debug,
     >     groups_maxx,points_maxx,points_maxx_group,particles_maxx,
     >     tweaks)
c
      implicit none
c
      integer groups_maxx,points_maxx,points_maxx_group,particles_maxx
      logical debug
      INTEGER m,n,lista(m),listb(m),nf(n)
      INTEGER j,k,l,tweaks
c
      do k=1,n
         nf(k)=k
      end do
c     
      do l=1,m
         j=lista(l)
         do while (nf(j) .ne. j)
            j=nf(j)
         end do
c     
         k=listb(l)
         do while (nf(k) .ne. k)
            k=nf(k)
         end do
c     
         if(j .ne. k) nf(j)=k
      end do
c     
      do j=1,n
         do while (nf(j) .ne. nf(nf(j)))
            nf(j)=nf(nf(j))
         end do
      end do
c
      if(debug) then
         do k=1,n
            write(18,*)k,nf(k)
         end do
      end if
      return
      end
