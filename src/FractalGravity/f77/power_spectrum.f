      subroutine power_spectrum(length,rho,lev,points_maxx)
c
      implicit none
c
      include 'maxx.inc'
c
      integer points_maxx
      integer wave,k,length,n,n_x,n_y,n_z,k_x,k_y,k_z
      integer k_nyq,k2,lev
      real rho(points_maxx),ak,ampl,scale
      real sum_0(1024),sum_1(1024),sum_2(1024)
c
      scale=2.0**lev
      k_nyq=length/2
      do k=1,length
         sum_2(k)=0.0
         sum_1(k)=0.0
         sum_0(k)=1.0e-30
      end do
c
      n=-1
      do n_z=1,length
         k_z=n_z-1
         if(k_z .gt. k_nyq)k_z=k_z-length
         do n_y=1,length
            k_y=n_y-1
            if(k_y .gt. k_nyq)k_y=k_y-length
            do n_x=1,length+2,2
               k_x=n_x/2
               k2=k_z**2+k_y**2+k_x**2
               ak=sqrt(float(k2))+1.0e-10
               wave=ak
               ampl=rho(n+2)**2+rho(n+3)**2
               if(wave .gt. 0 .and. ampl .gt. 0.0) then
                  n=n+2
                  sum_0(wave)=sum_0(wave)+1.0
                  sum_1(wave)=sum_1(wave)+ak
                  sum_2(wave)=sum_2(wave)+ampl
c                  write(17,*)'bb ',n,wave,ak,sum_0(wave),
c     >                 sum_1(wave),sum_2(wave)
               else
                  n=n+2
               end if
            end do
         end do
      end do
c
      do k=1,k_nyq
         write(18,19)k,scale*sum_1(k)/sum_0(k),
     >        sum_2(k)/sum_0(k),sum_0(k)
 19      format(i4,4(1pe11.3))
      end do
      return
      end
