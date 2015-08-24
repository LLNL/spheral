      program make_sphere
c     
      implicit none
c     
      integer n_tot,idum,n,points,p,maxx
      parameter (maxx=300000)
c     
      integer hlg(maxx)
      real x_0,y_0,z_0,r_0,alpha,sigma,expo,r_1,r_2,a
      real v_x,v_y,v_z,mass
      real mass_tot,pi,time,r,phi,c_phi,s_phi,c_theta,s_theta
      real step_length,dm,ran1,x_1,y_1,z_1,v_1,v_2,v_3
      real pos_x(maxx),pos_y(maxx),pos_z(maxx)
      real vel_x(maxx),vel_y(maxx),vel_z(maxx)
      real particle_mass(maxx),omega_1,omega_2,omega_3
c     
      print*,'x_0,y_0,z_0,r_0,alpha(<0),sigma'
      read(*,*)x_0,y_0,z_0,r_0,alpha,sigma
      print*,'omega_1,omega_2,omega_3'
      read(*,*)omega_1,omega_2,omega_3
      print*,'n_tot,mass_tot,step_length'
      read(*,*)n_tot,mass_tot,step_length
      print*,'idum'
      read(*,*)idum
c     
      pi=3.141592653589793
      expo=1.0/(3.0+alpha)
      time=0.0
c     
      dm=mass_tot/float(n_tot)
      do n=1,n_tot
        phi=2.0*pi*ran1(idum)
        c_phi=cos(phi)
        s_phi=sin(phi)
        c_theta=2.0*ran1(idum)-1.0
        s_theta=sqrt(1.0-c_theta**2)
        r=r_0*ran1(idum)**(1.0/(alpha+3.0))
        x_1=r*s_theta*c_phi
        pos_x(n)=x_1+x_0
        y_1=r*s_theta*s_phi
        pos_y(n)=y_1+y_0
        z_1=r*c_theta
        pos_z(n)=z_1+z_0
        v_1=z_1*omega_2-y_1*omega_3
        v_2=x_1*omega_3-z_1*omega_1
        v_3=y_1*omega_1-x_1*omega_2
c     
        r_1=ran1(idum)
        r_2=ran1(idum)
        vel_x(n)=sigma*sqrt(-2.0*alog(r_1))*cos(2.0*pi*r_2)+v_1
        vel_y(n)=sigma*sqrt(-2.0*alog(r_1))*sin(2.0*pi*r_2)+v_2
        r_1=ran1(idum)
        r_2=ran1(idum)
        vel_z(n)=sigma*sqrt(-2.0*alog(r_1))*cos(2.0*pi*r_2)+v_3
        particle_mass(n)=dm
      end do
c
      print*,'points,mass,x_0,y_0,z_0,v_x,v_y,v_z'
      read(*,*)points,mass,x_0,y_0,z_0,v_x,v_y,v_z
      if(points .gt. 0) then
        do p=1,points
          n_tot=n_tot+1
          pos_x(n_tot)=x_0
          pos_y(n_tot)=y_0
          pos_z(n_tot)=z_0
          vel_x(n_tot)=v_x
          vel_y(n_tot)=v_y
          vel_z(n_tot)=v_z
          particle_mass(n_tot)=mass
        end do
      end if
c     
      print*,'points,mass,x_0,y_0,z_0,v_x,v_y,v_z'
      read(*,*)points,mass,x_0,y_0,z_0,v_x,v_y,v_z
      if(points .gt. 0) then
        do p=1,points
          n_tot=n_tot+1
          pos_x(n_tot)=x_0
          pos_y(n_tot)=y_0
          pos_z(n_tot)=z_0
          vel_x(n_tot)=v_x
          vel_y(n_tot)=v_y
          vel_z(n_tot)=v_z
          particle_mass(n_tot)=mass
          hlg(n_tot)=1
        end do
      end if
c     
      a=0.0
      write(1)time,step_length,n_tot,a,a,a,1.0
      write(1)(pos_x(n),n=1,n_tot)
      write(1)(pos_y(n),n=1,n_tot)
      write(1)(pos_z(n),n=1,n_tot)
      write(1)(vel_x(n),n=1,n_tot)
      write(1)(vel_y(n),n=1,n_tot)
      write(1)(vel_z(n),n=1,n_tot)
      write(1)(particle_mass(n),n=1,n_tot)
      write(1)(hlg(n),n=1,n_tot)
      write(1)1,-1
c     
      stop 'it is all over'
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
