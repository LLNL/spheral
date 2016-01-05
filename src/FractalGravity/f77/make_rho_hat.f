      subroutine make_rho_hat(rho,level,length,points_maxx,
     >  parameters_integer,parameters_real,parameters_logical)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer points_maxx
      integer length,k_nyq,n,n_z,n_y,n_x
      integer k,k_z,k_y,k_x,k2,idum,length_1
      integer parameters_integer(1024),cut_wave
      integer spectrum_number,level_max_zel,level
c
      real ran1_tmp,amplitude,angle,power,scaling,q
      real rho(points_maxx),wave,ampl,ampl_0_no_cut
      real power_slope,cut_off,j_1,y,wave_0,pi,twopi,ampl0
      real filter,parameters_real(1024),vel_rho,step_wave
      real variance_0,variance_0_f,variance_2,variance_2_f
      real variance_0_no_cut,variance_2_no_cut
      real f_1(0:grid_length_max),f_2(0:grid_length_max)
c
      logical gauss_cut,parameters_logical(1024)
      logical short_x,short_y,short_z,short,shorty
      logical nyquist,real_x,real_y,real_z,really
c     
      nyquist(n_x,n_y)=n_x .eq. 0 .or. n_x .eq. n_y
      shorty(n_x,n_y)=abs(n_x) .le. n_y
      j_1(y)=3.0*(sin(y)-y*cos(y))/y**3
c     
      scaling=1.0
c
      pi=atan(1.0)*4.0
      twopi=atan(1.0)*8.0
c
      do k=0,length
        f_1(k)=(sin(2.0*pi*float(k)/float(length))*0.5)**2
        f_2(k)= (sin(pi*float(k)/float(length)))**2+1.0e-30
      end do
c
      power_slope=parameters_real(1)
      cut_off=parameters_real(2)
      wave_0=parameters_real(3)
      scaling=parameters_real(4)
      idum=-abs(parameters_integer(100+level))
c      idum=-abs(parameters_integer(1))
      spectrum_number=parameters_integer(2)
      level_max_zel=parameters_integer(3)
      gauss_cut=parameters_logical(1)
c     
      print*,'IDUM= ',idum
c
      step_wave=2.0**level
      if(level .eq. 0) then
         cut_wave=-1
      else
         cut_wave=length/4
      end if
c
      variance_0=0.0
      variance_2=0.0
      variance_0_no_cut=0.0
      variance_2_no_cut=0.0
      variance_0_f=0.0
      variance_2_f=0.0

c
      k_nyq=length/2
      length_1=length+1
      n=-1
c     
      do n_z=1,length
         k_z=n_z-1
         real_z=nyquist(k_z,k_nyq)
         if(k_z .gt. k_nyq)k_z=k_z-length
         short_z=shorty(k_z,cut_wave)
         do n_y=1,length
            k_y=n_y-1
            real_y=nyquist(k_y,k_nyq)
            if(k_y .gt. k_nyq)k_y=k_y-length
            short_y=shorty(k_y,cut_wave)
            do n_x=1,length+2,2
               k_x=n_x/2
               real_x=nyquist(k_x,k_nyq)
               really=real_x .and. real_y .and. real_z
               short_x=shorty(k_x,cut_wave)
               short=short_x .and. short_y .and. short_z
               k2=k_z**2+k_y**2+k_x**2
               wave=sqrt(float(k2)+1.0e-10)*step_wave
c     
               if(wave .gt. 0.001) then
                  q=wave/scaling
                  power=amplitude(q,power_slope,cut_off/scaling,
     >                 spectrum_number)
c                  write(17,*)wave,q,power
c      
                  ampl0=sqrt(power)
                  ampl=ampl0*sqrt(-2.0*alog(ran1_tmp(idum)))
                  ampl_0_no_cut=ampl0
               else
                  ampl=0.0
                  ampl0=0.0
                  ampl_0_no_cut=0.0
               end if
c
               if(short) then
                  ampl=0.0
                  ampl0=0.0
               end if
c
               if(k_y .eq. 0 .and. k_z .eq. 0) write(21,*)
     >              wave,ampl,ampl0
c     
               if(wave/wave_0 .gt. 0.01) then
                  if(gauss_cut) then
                     filter=exp(-0.5*(wave/wave_0)**2)
                  else
                     filter=j_1(wave/wave_0)
                  end if
               else
                  filter=1.0
               end if
c     
               vel_rho=(f_1(abs(k_x))+f_1(abs(k_y))+f_1(abs(k_z)))/
     >              (f_2(abs(k_x))+f_2(abs(k_y))+f_2(abs(k_z)))**2
c     
               if(really) then
                 angle=0.0
               else
                 angle=twopi*ran1_tmp(idum)
               end if
c
               n=n+2
c     
               rho(n)=ampl*cos(angle)
               rho(n+1)=ampl*sin(angle)
               variance_0=variance_0+ampl**2
               variance_0_no_cut=variance_0_no_cut+ampl_0_no_cut**2
               variance_0_f=variance_0_f+ampl**2*filter*filter
               variance_2=variance_2+ampl**2*vel_rho
               variance_2_f=variance_2_f+ampl**2*filter*filter*vel_rho
               variance_2_no_cut=variance_2_no_cut+ampl_0_no_cut**2*
     >              vel_rho
c     
c               write(22,22)n_x,n_y,n_z,k_x,k_y,k_z,n,rho(n),rho(n+1)
 22            format(6i4,i8,2(1pe13.5))
            end do
         end do
      end do
c

      if(level .eq. 0) then
         parameters_real(7)=sqrt(variance_0_f/variance_0)
         parameters_real(6)=sqrt(variance_2_f/variance_2)
         print*,'variances ',parameters_real(7),parameters_real(6)
      end if
c
      parameters_real(10+level*2)=variance_0_no_cut
      parameters_real(11+level*2)=variance_2_no_cut
c
      print*,'var_2 lev ',level,variance_2
c   
      return
      end
c
      real function amplitude(q,power_slope,cut_off,spectrum_number)
c
      implicit none
c
      real t,q,power_slope,cut_off
      integer spectrum_number
c
      amplitude=q**power_slope*exp(-(q/cut_off)**2)
c
      if(spectrum_number .eq. 0) then
         return
      else if(spectrum_number .eq. 1) then
         t=alog(1.0+2.34*q)/(2.34*q)*
     >        (1.0+13.0*q+(10.5*q)**2+(10.4*q)**3+(6.51*q)**4)**(-0.25)
         amplitude=amplitude*t*t
      end if
      return
      end
c
