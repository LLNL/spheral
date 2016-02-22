      subroutine list_buffer(point_0,up_x,up_y,up_z,
     >  down_x,down_y,down_z,corner,it_is_high,
     >  extra,
     >  points_maxx)
c     
      implicit none
c     
      integer points_maxx
c     
      integer point,point_0,corner,extra(7),add_points,p
      integer up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx)
      integer down_y(points_maxx),down_z(points_maxx)
c     
      logical it_is_high(-1:points_maxx)
c     
c     print*,'in list buffer'
      point=point_0
c     
      if(corner .eq. 1) then
         p=add_points(down_x,down_y,up_x,down_z,down_x,up_y,up_x,
     >        point,extra,points_maxx)
      else if(corner .eq. 2) then
         p=add_points(up_x,down_y,down_x,down_z,up_x,up_y,down_x,
     >        point,extra,points_maxx)
      else if(corner .eq. 3) then
         p=add_points(down_x,up_y,up_x,down_z,down_x,down_y,up_x,
     >        point,extra,points_maxx)
      else if(corner .eq. 4) then
        p=add_points(up_x,up_y,down_x,down_z,up_x,down_y,down_x,
     >        point,extra,points_maxx)
         else if(corner .eq. 5) then
         p=add_points(down_x,down_y,up_x,up_z,down_x,up_y,up_x,
     >        point,extra,points_maxx)
      else if(corner .eq. 6) then
         p=add_points(up_x,down_y,down_x,up_z,up_x,up_y,down_x,
     >        point,extra,points_maxx)
         else if(corner .eq. 7) then
         p=add_points(down_x,up_y,up_x,up_z,down_x,down_y,up_x,
     >        point,extra,points_maxx)
         else if(corner .eq. 8) then
        p=add_points(up_x,up_y,down_x,up_z,up_x,down_y,down_x,
     >        point,extra,points_maxx)
      end if
      return
      end
c
      integer function add_points(p1,p2,p3,p4,p5,p6,p7,p,extra,long)
c
      implicit none
c
      integer long,p
      integer p1(long),p2(long),p3(long),p4(long)
      integer p5(long),p6(long),p7(long),extra(7)
c
      p=p1(p)
      extra(1)=p
c
      p=p2(p)
      extra(2)=p
c
      p=p3(p)
      extra(3)=p
c
      p=p4(p)
      extra(4)=p
c
      p=p5(p)
      extra(5)=p
c
      p=p6(p)
      extra(6)=p
c
      p=p7(p)
      extra(7)=p
c
      add_points=p
c
      return
      end
