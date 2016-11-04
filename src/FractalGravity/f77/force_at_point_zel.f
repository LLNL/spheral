      subroutine force_at_point_zel(group,hoc_points,next_points,
     >   f_x,f_y,f_z,
     >   up_x,up_y,up_z,
     >   down_x,down_y,down_z,
     >   pointer,grid_length,level,highest_level_init,inside,
     >   debug,groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,tweaks
      logical inside(points_maxx),debug
c     
      real conv,f_x(points_maxx),f_y(points_maxx)
      real f_z(points_maxx)
c     
      integer point,hoc_points(groups_maxx),group,grid_length
      integer level(groups_maxx),highest_level_init
      integer up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx)
      integer p_up_x,p_down_x,p_up_y,p_down_y,p_up_z,p_down_z,p_p
      integer down_x(points_maxx),down_y(points_maxx)
      integer down_z(points_maxx)
      integer pointer(points_maxx),next_points(points_maxx)
      integer p_p_up,p_p_down,p_p_up_up,p_p_down_down
      integer p_p_down_up,p_p_up_down
      integer p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8
      integer p_up_y_up_x,p_up_y_down_x,p_down_y_down_x,p_down_y_up_x
c     
      conv=float(grid_length)*2.0**(level(group)-1)
      if(debug)
     >   write(34,*)'enter force at point zel ',group,level(group),conv
c     
      point=hoc_points(group)
      do while(point .gt. 0)
         p_p=-pointer(point)
         p_up_x=up_x(point)
         p_up_y=up_y(point)
         p_up_z=up_z(point)
         p_down_x=down_x(point)
         p_down_y=down_y(point)
         p_down_z=down_z(point)
c     
c     if(inside(point)) then
c     
c         write(47,48)group,point,f_x(point),f_y(point),f_z(point)
 48      format(i4,i8,3(1pe13.5))
         if(p_p .lt. 0) then
            f_x(point)=f_x(point)+f_x(-p_p)
            f_y(point)=f_y(point)+f_y(-p_p)
            f_z(point)=f_z(point)+f_z(-p_p)
c     
         else if(p_p .eq. 2 .or. p_p .eq. 8 .or.
     >           p_p .eq. 20 .or. p_p .eq. 26) then
            p_p_up=pointer(p_up_x)
            p_p_down=pointer(p_down_x)
            f_x(point)=f_x(point)+0.5*(
     >           f_x(p_p_up)+f_x(p_p_down))
            f_y(point)=f_y(point)+0.5*(
     >           f_y(p_p_up)+f_y(p_p_down))
            f_z(point)=f_z(point)+0.5*(
     >           f_z(p_p_up)+f_z(p_p_down))
c     
         else if(p_p .eq. 4 .or. p_p .eq. 6 .or. 
     >           p_p .eq. 22 .or. p_p .eq. 24) then
            p_p_up=pointer(p_up_y)
            p_p_down=pointer(p_down_y)
            f_x(point)=f_x(point)+0.5*(
     >           f_x(p_p_up)+f_x(p_p_down))
            f_y(point)=f_y(point)+0.5*(
     >           f_y(p_p_up)+f_y(p_p_down))
            f_z(point)=f_z(point)+0.5*(
     >           f_z(p_p_up)+f_z(p_p_down))
c     
         else if(p_p .eq. 10 .or. p_p .eq. 12 .or. 
     >           p_p .eq. 16 .or. p_p .eq. 18) then
            p_p_up=pointer(p_up_z)
            p_p_down=pointer(p_down_z)
            f_x(point)=f_x(point)+0.5*(
     >           f_x(p_p_up)+f_x(p_p_down))
            f_y(point)=f_y(point)+0.5*(
     >           f_y(p_p_up)+f_y(p_p_down))
            f_z(point)=f_z(point)+0.5*(
     >           f_z(p_p_up)+f_z(p_p_down))
c     
         else if(p_p .eq. 5 .or. p_p .eq. 23) then
            p_p_up_up=pointer(up_x(p_up_y))
            p_p_down_up=pointer(down_x(p_up_y))
            p_p_down_down=pointer(down_x(p_down_y))
            p_p_up_down=pointer(up_x(p_down_y))
            f_x(point)=f_x(point)+0.25*(
     >           f_x(p_p_down_down)+f_x(p_p_up_down)+
     >           f_x(p_p_down_up)+f_x(p_p_up_up))
            f_y(point)=f_y(point)+0.25*(
     >           f_y(p_p_down_down)+f_y(p_p_up_down)+
     >           f_y(p_p_down_up)+f_y(p_p_up_up))
            f_z(point)=f_z(point)+0.25*(
     >           f_z(p_p_down_down)+f_z(p_p_up_down)+
     >           f_z(p_p_down_up)+f_z(p_p_up_up))
         else if(p_p .eq. 11 .or. p_p .eq. 17) then
            p_p_up_up=pointer(up_x(p_up_z))
            p_p_down_up=pointer(down_x(p_up_z))
            p_p_down_down=pointer(down_x(p_down_z))
            p_p_up_down=pointer(up_x(p_down_z))
            f_x(point)=f_x(point)+0.25*(
     >           f_x(p_p_down_down)+f_x(p_p_up_down)+
     >           f_x(p_p_down_up)+f_x(p_p_up_up))
            f_y(point)=f_y(point)+0.25*(
     >           f_y(p_p_down_down)+f_y(p_p_up_down)+
     >           f_y(p_p_down_up)+f_y(p_p_up_up))
            f_z(point)=f_z(point)+0.25*(
     >           f_z(p_p_down_down)+f_z(p_p_up_down)+
     >           f_z(p_p_down_up)+f_z(p_p_up_up))
         else if(p_p .eq. 13 .or. p_p .eq. 15) then
            p_p_up_up=pointer(up_y(p_up_z))
            p_p_down_up=pointer(down_y(p_up_z))
            p_p_down_down=pointer(down_y(p_down_z))
            p_p_up_down=pointer(up_y(p_down_z))
            f_x(point)=f_x(point)+0.25*(
     >           f_x(p_p_down_down)+f_x(p_p_up_down)+
     >           f_x(p_p_down_up)+f_x(p_p_up_up))
            f_y(point)=f_y(point)+0.25*(
     >           f_y(p_p_down_down)+f_y(p_p_up_down)+
     >           f_y(p_p_down_up)+f_y(p_p_up_up))
            f_z(point)=f_z(point)+0.25*(
     >           f_z(p_p_down_down)+f_z(p_p_up_down)+
     >           f_z(p_p_down_up)+f_z(p_p_up_up))
         else if(p_p .eq. 14) then
            p_up_y_up_x=up_y(p_up_x)
            p_down_y_up_x=down_y(p_up_x)
            p_up_y_down_x=up_y(p_down_x)
            p_down_y_down_x=down_y(p_down_x)
            p_2=pointer(down_z(p_down_y_up_x))
            p_6=pointer(up_z(p_down_y_up_x))
            p_4=pointer(down_z(p_up_y_up_x))
            p_8=pointer(up_z(p_up_y_up_x))
            p_1=pointer(down_z(p_down_y_down_x))
            p_5=pointer(up_z(p_down_y_down_x))
            p_3=pointer(down_z(p_up_y_down_x))
            p_7=pointer(up_z(p_up_y_down_x))
c     
            f_x(point)=f_x(point)+0.125*(
     >           f_x(p_1)+f_x(p_2)+f_x(p_3)+f_x(p_4)+
     >           f_x(p_5)+f_x(p_6)+f_x(p_7)+f_x(p_8))
            f_y(point)=f_y(point)+0.125*(
     >           f_y(p_1)+f_y(p_2)+f_y(p_3)+f_y(p_4)+
     >           f_y(p_5)+f_y(p_6)+f_y(p_7)+f_y(p_8))
            f_z(point)=f_z(point)+0.125*(
     >           f_z(p_1)+f_z(p_2)+f_z(p_3)+f_z(p_4)+
     >           f_z(p_5)+f_z(p_6)+f_z(p_7)+f_z(p_8))
         else 
            stop 'force not with you at this point'
         end if
c         write(57,48)group,point,f_x(point),f_y(point),f_z(point)
c     end if
         point=next_points(point)
      end do
c     
      return
      end
