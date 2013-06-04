      subroutine find_neighbors(point,adj,logic,
     >  up_x,up_y,up_z,down_x,down_y,down_z,points_maxx)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer points_maxx
      integer point,i,j,k,up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx)
      integer down_y(points_maxx),down_z(points_maxx)
      integer adj(-1:1,-1:1,-1:1),g_1,g_2,g_3
      logical logic(-1:points_maxx)
c     
      do k=-1,1
        do j=-1,1
          do i=-1,1
            adj(i,j,k)=-1
          end do
        end do
      end do
c     
      adj(0,0,0)=point
c     
      adj(1,0,0)=g_1(point,up_x,logic)
      adj(0,1,0)=g_1(point,up_y,logic)
      adj(0,0,1)=g_1(point,up_z,logic)
      adj(-1,0,0)=g_1(point,down_x,logic)
      adj(0,-1,0)=g_1(point,down_y,logic)
      adj(0,0,-1)=g_1(point,down_z,logic)
c
      adj(-1,-1,0)=g_2(point,down_x,down_y,logic)
      adj(1,-1,0)=g_2(point,up_x,down_y,logic)
      adj(-1,1,0)=g_2(point,down_x,up_y,logic)
      adj(1,1,0)=g_2(point,up_x,up_y,logic)
      adj(1,0,1)=g_2(point,up_x,up_z,logic)
      adj(0,1,1)=g_2(point,up_y,up_z,logic)
      adj(-1,0,1)=g_2(point,down_x,up_z,logic)
      adj(0,-1,1)=g_2(point,down_y,up_z,logic)
      adj(1,0,-1)=g_2(point,up_x,down_z,logic)
      adj(0,1,-1)=g_2(point,up_y,down_z,logic)
      adj(-1,0,-1)=g_2(point,down_x,down_z,logic)
      adj(0,-1,-1)=g_2(point,down_y,down_z,logic)
c     
      adj(1,1,1)=g_3(point,up_x,up_y,up_z,logic)
      adj(-1,1,1)=g_3(point,down_x,up_y,up_z,logic)
      adj(1,-1,1)=g_3(point,up_x,down_y,up_z,logic)
      adj(-1,-1,1)=g_3(point,down_x,down_y,up_z,logic)
      adj(1,1,-1)=g_3(point,up_x,up_y,down_z,logic)
      adj(-1,1,-1)=g_3(point,down_x,up_y,down_z,logic)
      adj(1,-1,-1)=g_3(point,up_x,down_y,down_z,logic)
      adj(-1,-1,-1)=g_3(point,down_x,down_y,down_z,logic)
c     
      return
      end
c
      integer function g_1(i,dir_1,logic)
c
      implicit none
c
      integer i,dir_1(*)
      logical logic(-1:*)
c
      g_1=-1
      if(i .lt. 0) return
      g_1=dir_1(i)
      if(.not. logic(g_1)) g_1=-1
c
      return
      end
c
      integer function g_2(i,dir_1,dir_2,logic)
c
      integer i,dir_1(*),dir_2(*),g_1
      logical logic(-1:*)
c
       g_2=g_1(g_1(i,dir_1,logic),dir_2,logic)
      if(g_2 .gt. 0) return
      g_2=g_1(g_1(i,dir_2,logic),dir_1,logic)
c
      return
      end
c
      integer function g_3(i,dir_1,dir_2,dir_3,logic)
c
      integer i,dir_1(*),dir_2(*),dir_3(*),g_2,g_1
      logical logic(-1:*)
c
      g_3=g_1(g_2(i,dir_1,dir_2,logic),dir_3,logic)
      if(g_3 .gt. 0) return
      g_3=g_1(g_2(i,dir_2,dir_3,logic),dir_1,logic)
      if(g_3 .gt. 0) return
      g_3=g_1(g_2(i,dir_1,dir_3,logic),dir_2,logic)
c
      return
      end
