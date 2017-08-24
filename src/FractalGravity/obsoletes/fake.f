      program fake
c     
      implicit none
c     
      logical periodic,debug
c     
      real pos_x(50000),pos_y(50000),pos_z(50000)
      real x0,y0,z0,r0,r,angle,twopi,ran1,cos_ang,sin_ang
      real force_x(50000),force_y(50000),force_z(50000)
      real potential(50000),particle_mass(50000),dm
c     
      integer number_particles,grid_length,minimum_number,level_max
      integer moat,idum,total,block,blocks,part,particles
      integer random_offset,tweaks,memory_value
c     
      print*,'periodic,debug t/f,idum ,offset(I),tweaks(I)'
      read(5,*)periodic,debug,idum,random_offset,tweaks
      print*,'minimum number,grid_length,moat,level_max'
      read(5,*)minimum_number,grid_length,moat,level_max
c     
      twopi=8.0*atan(1.0)
      total=0
      print*,'number of blocks'
      read(5,*)blocks
      do block=1,blocks
       print*,block,total
       print*,'particles,x0,y0,z0,r0,dm'
       read(5,*)particles,x0,y0,z0,r0,dm
       do part=1,particles
        total=total+1
        r=r0*ran1(idum)
        angle=twopi*ran1(idum)
        cos_ang=(ran1(idum)-0.5)*1.999
        sin_ang=sqrt(1.0-cos_ang**2)
        pos_x(total)=x0+r*cos(angle)*sin_ang
        pos_y(total)=y0+r*sin(angle)*sin_ang
        pos_z(total)=z0+r*cos_ang
        particle_mass(total)=dm
       end do
      end do
      number_particles=total
c     
      call fractal_gravity(number_particles,grid_length,periodic,
     >  minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >  potential,force_x,force_y,force_z,
     >  moat,random_offset,debug,tweaks,memory_value)
c     
 10   format(i8,2f10.4)
      stop 'it is all over fake'
      end
c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
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
      ran1=min(AM*iy,RNMX)
      return
      END
