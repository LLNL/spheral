      program pancake
c     
      implicit none
c     
      integer number_particles,idum,n,n_t,n_1,n_x,n_y,n_z,n_pan,n_p
      real dt,omega_0,omega_lambda,delta
      real k_1(10),k_2(10),k_3(10),k2,a_cr(10),pi,dm
      real x_0,y_0,z_0,coef,arg,time,x_c(10),y_c(10),z_c(10),
     >  ran1,arad,offset
      real particle_mass(300000),pos_x(300000)
      real pos_y(300000),pos_z(300000)
      real vel_x(300000),vel_y(300000),vel_z(300000)
      real d_x,d_y,d_z,pexp,fact
      print*,'number_particles,dt,omega_0,omega_lambda,idum'
      read(*,*)number_particles,dt,omega_0,omega_lambda,idum
c     
      print*,'number of pancakes'
      read(*,*)n_pan
      do n=1,n_pan
       print*,'k_1,k_2,k_3,a_cr'
       read(*,*)k_1(n),k_2(n),k_3(n),a_cr(n)
       print*,'x_c,y_c,z_c'
       read(*,*)x_c(n),y_c(n),z_c(n)
      end do
c     
      pi=3.1415926535
      pexp=1.0
      fact=1.0/sqrt(1.0-dt)
c     
      n_1=(float(number_particles)+0.1)**(1.0/3.0)
      delta=1.0/float(n_1)
      offset=delta/2.0
      if(idum .eq. 1) number_particles=n_1**3
      dm=3.0/(8.0*pi)*omega_0/float(number_particles)
      n_t=0
      do n=1,number_particles
       if(idum .ne. 0) then
        x_0=ran1(idum)
        y_0=ran1(idum)
        z_0=ran1(idum)
       else
        n_z=n_t/(n_1*n_1)
        n_y=mod(n_t/n_1,n_1)
        n_x=mod(n_t,n_1)
        x_0=float(n_x)*delta+offset
        y_0=float(n_y)*delta+offset
        z_0=float(n_z)*delta+offset
        n_t=n_t+1
       end if
c     
       particle_mass(n)=dm
c     
       d_x=0.0
       d_y=0.0
       d_z=0.0
       do n_p=1,n_pan
        k2=k_1(n_p)**2+k_2(n_p)**2+k_3(n_p)**2
        coef=1.0/(2.0*pi*k2*a_cr(n_p))
        arg=sin(2.0*pi*(k_1(n_p)*(x_0-x_c(n_p))+
     >    k_2(n_p)*(y_0-y_c(n_p))+
     >    k_3(n_p)*(z_0-z_c(n_p))))
        d_x=d_x-k_1(n_p)*coef*arg
        d_y=d_y-k_2(n_p)*coef*arg
        d_z=d_z-k_3(n_p)*coef*arg
       end do
       pos_x(n)=x_0+d_x
       pos_y(n)=y_0+d_y
       pos_z(n)=z_0+d_z
       vel_x(n)=d_x*fact
       vel_y(n)=d_y*fact
       vel_z(n)=d_z*fact
      end do
c     
      time=2.0/3.0
      arad=1.0
      write(2)time,dt,number_particles,
     >  arad,omega_0,omega_lambda,pexp
      write(2)(pos_x(n),n=1,number_particles)
      write(2)(pos_y(n),n=1,number_particles)
      write(2)(pos_z(n),n=1,number_particles)
      write(2)(vel_x(n),n=1,number_particles)
      write(2)(vel_y(n),n=1,number_particles)
      write(2)(vel_z(n),n=1,number_particles)
      write(2)(particle_mass(n),n=1,number_particles)
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
