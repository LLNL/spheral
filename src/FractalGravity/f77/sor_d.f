      subroutine sor_d(group,hoc_points,next_points,
     >     n_pointer,pointer_left_x,
     >     pot,rho,up_x,up_y,up_z,
     >     down_x,down_y,down_z,inside,grav_const,
     >     rjac,eps,memory_value,debug,maxits,
     >     groups_maxx,points_maxx)
c     
      implicit none
c     
      include 'maxes.inc'
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx,point,n
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer pointer_left_x(8*grid_length_max**2)
      integer up_x(points_maxx),up_y(points_maxx),up_z(points_maxx)
      integer down_x(points_maxx),down_y(points_maxx)
      integer down_z(points_maxx),point_left,ipass,maxits,n_points
      integer maxits_real,n_inside,memory_value
      integer n_pointer,pointer(points_max)
c     
      logical inside(points_maxx),debug,it_is_inside
c     
      real pot(points_maxx),rho(points_maxx)
      real grav_const,rjac,eps
c     
      double precision pot_d(points_max_group),omega_6,anorm,anorm2
      double precision rho_d(points_max_group),resid,anormf,anorm2f
      double precision eps_d,rjac_d,omega
c     
      eps_d=eps
      rjac_d=rjac
      anormf=0.0
      n_points=0
      n_inside=0
      point=hoc_points(group)
      do while (point .gt. 0)
         n_points=n_points+1
         pot_d(n_points)=pot(point)
         rho_d(n_points)=rho(point)*grav_const
         if(inside(point)) then
            anormf=anormf+abs(pot_d(n_points))
            n_inside=n_inside+1
         end if
         pointer(point)=n_points
         point=next_points(point)
      end do
c     
      maxits_real=maxits
      if(n_inside .eq. 1) maxits_real=1
c     
      anorm2f=3.0*anormf/float(n_inside)
c     
      omega=1.0
      omega_6=omega/6.0
      do n=1,maxits_real
         anorm=0.0
         anorm2=-1.0
         do ipass=1,2
            do point_left=n_pointer,1,-1
               point=pointer_left_x(point_left)
               if(ipass .eq. 2)point=up_x(point)
               do while(it_is_inside(point,inside))
                  resid=
     >                 pot_d(pointer(up_x(point)))+
     >                 pot_d(pointer(up_y(point)))+
     >                 pot_d(pointer(up_z(point)))+
     >                 pot_d(pointer(down_x(point)))+
     >                 pot_d(pointer(down_y(point)))+
     >                 pot_d(pointer(down_z(point)))-
     >                 6.0*pot_d(pointer(point))-rho_d(pointer(point))
                  anorm=anorm+abs(resid)
                  anorm2=max(anorm2,abs(resid))
                  pot_d(pointer(point))=pot_d(pointer(point))+
     >                 omega_6*resid
c     
                  point=up_x(point)
                  if(point .gt. 0) point=up_x(point)
               end do
            end do
            if(n .eq. 1 .and. ipass .eq. 1) then
               omega=1.0/(1.0-0.5*rjac_d**2)
            else
               omega=1.0/(1.0-0.25*rjac_d**2*omega)
            end if
            omega_6=omega/6.0
         end do
         if(anorm .lt. eps_d*anormf) then
            if(anorm2 .lt. eps_d*anorm2f) then
               print 1,group,n,n_inside,
     >              anorm,anormf,anorm2,anorm2f,rjac
c     
               point=hoc_points(group)
               do while (point .gt. 0)
                  if(inside(point)) pot(point)=pot_d(pointer(point))
                  point=next_points(point)
               end do
               return
            end if
         end if
      end do
c
      if(maxits_real .eq. 1) then
         print 1,group,n,n_inside,
     >        anorm,anormf,anorm2,anorm2f,rjac
         return
      else
         print 2,group,n,n_inside,
     >        anorm,anormf,anorm2,anorm2f,rjac
         memory_value=8
         print*,'memory value= ',memory_value
         return
      end if
c     
 1    format('sor conv ',3i8,5(1pe12.4))
 2    format('sor not conv ',3i8,5(1pe12.4))
      stop
      end
