      subroutine offset(what,random_offset,number_particles,
     >  grid_multiply,pos_x,pos_y,pos_z,tweaks)
c     
      implicit none
c     
      integer what,number_particles,particle,grid_multiply
      integer random_offset,idum,tweaks
      real pos_x(number_particles),pos_y(number_particles)
      real pos_z(number_particles)
      real offset_x,offset_y,offset_z,ran1_tmp
      data offset_x,offset_y,offset_z,idum/3*0.0,0/
c     
c***  offset is run to avoid badness if the data is on a uniform
c***  lattice coincident with the grid points on level 0.
c***  the points are offset one way at the start of the program 
c***  and are moved back afterwards 
c     
c     
      if(what .gt. 0) then
        if(random_offset .eq. 0) then
          offset_x=0.0
          offset_y=0.0
          offset_z=0.0
          return
        else
          idum=-abs(random_offset)
          offset_x=(2.0*(ran1_tmp(idum)-0.5))/float(grid_multiply)
          offset_y=(2.0*(ran1_tmp(idum)-0.5))/float(grid_multiply)
          offset_z=(2.0*(ran1_tmp(idum)-0.5))/float(grid_multiply)
        end if
c     
        do particle=1,number_particles
          pos_x(particle)=pos_x(particle)+offset_x
          pos_y(particle)=pos_y(particle)+offset_y
          pos_z(particle)=pos_z(particle)+offset_z
        end do
c     
      else if(what .lt. 0) then
c     
        do particle=1,number_particles
          pos_x(particle)=pos_x(particle)-offset_x
          pos_y(particle)=pos_y(particle)-offset_y
          pos_z(particle)=pos_z(particle)-offset_z
        end do
      else
        stop 'die in offset'
      end if
c     
      print*,' offset= ',offset_x,offset_y,offset_z
c     
      return
      end
c
      FUNCTION ran1_tmp(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1_tmp,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1_tmp=min(AM*iy,RNMX)
      return
      END
