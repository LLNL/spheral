      subroutine sor(group,hoc_a,update_a,hoc_b,update_b,
     >  pot,rho,up_x,up_y,up_z,
     >  down_x,down_y,down_z,grav_const,
     >  rjac,eps,debug,maxits,points_maxx)
c     
      implicit none
c     
      integer points_maxx,point,n
      integer group,maxits_real
      integer hoc_a(2),update_a(points_maxx)
      integer hoc_b(2),update_b(points_maxx)
      integer up_x(points_maxx),up_y(points_maxx),up_z(points_maxx)
      integer down_x(points_maxx),down_y(points_maxx)
      integer down_z(points_maxx),ipass,maxits,n_points
c     
      logical debug,first_time
c     
      real pot(points_maxx),rho(points_maxx)
      real grav_const,rjac,eps,sixi,anormf,anorm2f
      real omega,anorm,resid,anorm2
c     
      sixi=1.0/6.0
      anormf=0.0
      n_points=0
      do ipass=1,2
        point=hoc_a(ipass)
        do while(point .gt. 0)
          n_points=n_points+1
          anormf=anormf+abs(pot(point))
          point=update_a(point)
        end do
      end do
c     
      anorm2f=3.0*anormf/float(n_points)
c     
      maxits_real=maxits
      if(n_points .eq. 1) maxits_real=1
c     
      first_time=.true.
      omega=1.0
      do n=1,maxits_real
        anorm=0.0
        anorm2=-1.0
        do ipass=1,2
          point=hoc_a(ipass)
          do while(point .gt. 0)
            resid=
     >        pot(up_x(point))+pot(up_y(point))+
     >        pot(up_z(point))+pot(down_x(point))+
     >        pot(down_y(point))+pot(down_z(point))-
     >        6.0*pot(point)-grav_const*rho(point)
            anorm=anorm+abs(resid)
            anorm2=max(anorm2,abs(resid))
            pot(point)=pot(point)+omega*resid/6.0
            point=update_a(point)
          end do
          if(first_time) then
            omega=1.0/(1.0-0.5*rjac**2)
            first_time=.false.
          else
            omega=1.0/(1.0-0.25*rjac**2*omega)
          end if
        end do
        if(n_points .gt. 50000)print*,'sor test ',
     >    n,anorm,anormf,anorm2,anorm2f
        if(anorm .lt. eps*anormf) then
          if(anorm2 .lt. eps*anorm2f) then
            print 1,group,n,n_points,anorm,anormf,anorm2,anorm2f
            return
          end if
        end if
c     
        if(maxits_real .eq. 1) return
c     
        anorm=0.0
        anorm2=-1.0
        do ipass=1,2
          point=hoc_b(ipass)
          do while(point .gt. 0)
            resid=
     >        pot(up_x(point))+pot(up_y(point))+
     >        pot(up_z(point))+pot(down_x(point))+
     >        pot(down_y(point))+pot(down_z(point))-
     >        6.0*pot(point)-grav_const*rho(point)
            anorm=anorm+abs(resid)
            anorm2=max(anorm2,abs(resid))
            pot(point)=pot(point)+omega*resid/6.0
            point=update_b(point)
          end do
          omega=1.0/(1.0-0.25*rjac**2*omega)
        end do
        if(n_points .gt. 50000)print*,'sor test ',
     >    n,anorm,anormf,anorm2,anorm2f
        if(anorm .lt. eps*anormf) then
          if(anorm2 .lt. eps*anorm2f) then
            print 1,group,n,n_points,anorm,anormf,anorm2,anorm2f
            return
          end if
        end if
      end do
c     
 1    format('sor conv ',3i8,4(1pe12.4))
      print*,'not converged ',group,anorm,anormf,n,n_points
      stop 'not converged'
      end
