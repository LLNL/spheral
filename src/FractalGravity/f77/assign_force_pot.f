      subroutine assign_force_pot(pot,highest_used_point,
     >     hoc_points,next_points,
     >     f_x,f_y,f_z,
     >     x,y,z,lev,highest_level_init,level_max,
     >     length,points_maxx,tweaks,
     >     parameters_integer,parameters_real,parameters_logical)
c
      implicit none
c
      integer points_maxx
      integer highest_used_point,lev,level_max,length
      integer parameters_integer(1024)
      integer x(points_maxx),y(points_maxx)
      integer z(points_maxx),next_points(points_maxx)
      integer hoc_points,p,p_x,p_y,p_z,point,wrapping
      integer division,n,p_x_up,p_y_up,p_z_up,tweaks
      integer p_x_down,p_y_down,p_z_down,highest_level_init
      integer n_x_up,n_x_down,n_y_up,n_y_down,n_z_up,n_z_down
c
      logical dont_do_it,parameters_logical(1024)
c
      real pot(points_maxx),conv,parameters_real(1024),amplify
      real f_x(points_maxx),f_y(points_maxx),f_z(points_maxx)
      real variance_2_no_cut_0,variance_2_no_cut
c
      n(p_x,p_y,p_z,length)=1+p_x+length*(p_y+p_z*length)
c
      wrapping=length*2**(level_max-lev)
      division=2**(level_max-lev)
      conv=float(length/2)
      variance_2_no_cut_0=parameters_real(11)
      variance_2_no_cut=parameters_real(11+lev*2)
      amplify=sqrt((variance_2_no_cut_0/variance_2_no_cut)*2.0**lev)
      write(68,*)'assign ',lev,parameters_real(11),
     >     parameters_real(11+2*lev),amplify
c
c     write(38,*)'assign pot ',lev,wrapping,division,
c    >     highest_used_point
c
      dont_do_it=lev .gt. highest_level_init
c
      point=hoc_points
      do while(point .gt. 0)
         if(dont_do_it) then
            pot(point)=0.0
            f_x(point)=0.0
            f_y(point)=0.0
            f_z(point)=0.0
         else
            p_x=(mod(x(point),wrapping))/division
            p_y=(mod(y(point),wrapping))/division
            p_z=(mod(z(point),wrapping))/division
            p=1+p_x+length*(p_y+p_z*length)
c     
            p_x_up=mod(p_x+length+1,length)
            p_y_up=mod(p_y+length+1,length)
            p_z_up=mod(p_z+length+1,length)
            p_x_down=mod(p_x+length-1,length)
            p_y_down=mod(p_y+length-1,length)
            p_z_down=mod(p_z+length-1,length)
c     
            n_x_up=n(p_x_up,p_y,p_z,length)+highest_used_point
            n_x_down=n(p_x_down,p_y,p_z,length)+highest_used_point
            n_y_up=n(p_x,p_y_up,p_z,length)+highest_used_point
            n_y_down=n(p_x,p_y_down,p_z,length)+highest_used_point
            n_z_up=n(p_x,p_y,p_z_up,length)+highest_used_point
            n_z_down=n(p_x,p_y,p_z_down,length)+highest_used_point
c     
            pot(point)=pot(p+highest_used_point)
            f_x(point)=(pot(n_x_down)-pot(n_x_up))*conv*amplify
            f_y(point)=(pot(n_y_down)-pot(n_y_up))*conv*amplify
            f_z(point)=(pot(n_z_down)-pot(n_z_up))*conv*amplify
         end if
c         write(38,*)'pot ',point,pot(point),f_x(point),
c     >        f_y(point),f_z(point)
         point=next_points(point)
      end do
c
      return
      end
