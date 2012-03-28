      subroutine find_neighbors(point,adj,
     >  up_x,up_y,up_z,down_x,down_y,down_z,points_maxx)
c     
      implicit none
c     
      integer points_maxx
      integer point,i,j,k,up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx)
      integer down_y(points_maxx),down_z(points_maxx)
      integer adj(-1:1,-1:1,-1:1),g_2,g_3,a
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
      adj(1,0,0)=up_x(point)
      adj(0,1,0)=up_y(point)
      adj(0,0,1)=up_z(point)
c     
      adj(-1,0,0)=down_x(point)
      adj(0,-1,0)=down_y(point)
      adj(0,0,-1)=down_z(point)
c     
      adj(-1,-1,0)=g_2(point,down_x,down_y)
      if(adj(-1,-1,0) .lt. 0) adj(-1,-1,0)=g_2(point,down_y,down_x)
c     
      adj(1,-1,0)=g_2(point,up_x,down_y)
      if(adj(1,-1,0) .lt. 0) adj(1,-1,0)=g_2(point,down_y,up_x)
c     
      adj(-1,1,0)=g_2(point,down_x,up_y)
      if(adj(-1,1,0) .lt. 0)adj(-1,1,0)=g_2(point,up_y,down_x)
c     
      adj(1,1,0)=g_2(point,up_x,up_y)
      if(adj(1,1,0) .lt. 0) adj(1,1,0)=g_2(point,up_y,up_x)
c     
      adj(1,0,1)=g_2(point,up_x,up_z)
      if(adj(1,0,1) .lt. 0) adj(1,0,1)=g_2(point,up_z,up_x)
c     
      adj(0,1,1)=g_2(point,up_y,up_z)
      if(adj(0,1,1) .lt. 0) adj(0,1,1)=g_2(point,up_z,up_y)
c     
      adj(-1,0,1)=g_2(point,down_x,up_z)
      if(adj(-1,0,1) .lt. 0) adj(-1,0,1)=g_2(point,up_z,down_x)
c     
      adj(0,-1,1)=g_2(point,down_y,up_z)
      if(adj(0,-1,1) .lt. 0) adj(0,-1,1)=g_2(point,up_z,down_y)
c     
      adj(1,0,-1)=g_2(point,up_x,down_z)
      if(adj(1,0,-1) .lt. 0) adj(1,0,-1)=g_2(point,down_z,up_x)
c     
      adj(0,1,-1)=g_2(point,up_y,down_z)
      if(adj(0,1,-1) .lt. 0) adj(0,1,-1)=g_2(point,down_z,up_y)
c     
      adj(-1,0,-1)=g_2(point,down_x,down_z)
      if(adj(-1,0,-1) .lt. 0) adj(-1,0,-1)=g_2(point,down_z,down_x)
c     
      adj(0,-1,-1)=g_2(point,down_y,down_z)
      if(adj(0,-1,-1) .lt. 0) adj(0,-1,-1)=g_2(point,down_z,down_y)
c     
      a=-1
      a=g_3(point,up_x,up_y,up_z)
      if(a .lt. 0) then
        a=g_3(point,up_x,up_z,up_y)
        if(a .lt. 0) then
          a=g_3(point,up_y,up_z,up_x)
          if(a .lt. 0) then
            a=g_3(point,up_y,up_x,up_z)
            if(a .lt. 0) then
              a=g_3(point,up_z,up_x,up_y)
              if(a .lt. 0) then
                a=g_3(point,up_z,up_y,up_x)
              end if
            end if
          end if
        end if
      end if
      adj(1,1,1)=a
c
      a=-1
      a=g_3(point,down_x,up_y,up_z)
      if(a .lt. 0) then
        a=g_3(point,down_x,up_z,up_y)
        if(a .lt. 0) then
          a=g_3(point,up_y,up_z,down_x)
          if(a .lt. 0) then
            a=g_3(point,up_y,down_x,up_z)
            if(a .lt. 0) then
              a=g_3(point,up_z,down_x,up_y)
              if(a .lt. 0) then
                a=g_3(point,up_z,up_y,down_x)
              end if
            end if
          end if
        end if
      end if
      adj(-1,1,1)=a
c
      a=-1
      a=g_3(point,up_x,down_y,up_z)
      if(a .lt. 0) then
        a=g_3(point,up_x,up_z,down_y)
        if(a .lt. 0) then
          a=g_3(point,down_y,up_z,up_x)
          if(a .lt. 0) then
            a=g_3(point,down_y,up_x,up_z)
            if(a .lt. 0) then
              a=g_3(point,up_z,up_x,down_y)
              if(a .lt. 0) then
                a=g_3(point,up_z,down_y,up_x)
              end if
            end if
          end if
        end if
      end if
      adj(1,-1,1)=a
c
      a=-1
      a=g_3(point,down_x,down_y,up_z)
      if(a .lt. 0) then
        a=g_3(point,down_x,up_z,down_y)
        if(a .lt. 0) then
          a=g_3(point,down_y,up_z,down_x)
          if(a .lt. 0) then
            a=g_3(point,down_y,down_x,up_z)
            if(a .lt. 0) then
              a=g_3(point,up_z,down_x,down_y)
              if(a .lt. 0) then
                a=g_3(point,up_z,down_y,down_x)
              end if
            end if
          end if
        end if
      end if
      adj(-1,-1,1)=a
c
      a=-1
      a=g_3(point,up_x,up_y,down_z)
      if(a .lt. 0) then
        a=g_3(point,up_x,down_z,up_y)
        if(a .lt. 0) then
          a=g_3(point,up_y,down_z,up_x)
          if(a .lt. 0) then
            a=g_3(point,up_y,up_x,down_z)
            if(a .lt. 0) then
              a=g_3(point,down_z,up_x,up_y)
              if(a .lt. 0) then
                a=g_3(point,down_z,up_y,up_x)
              end if
            end if
          end if
        end if
      end if
      adj(1,1,-1)=a
c     
      a=-1
      a=g_3(point,down_x,up_y,down_z)
      if(a .lt. 0) then
        a=g_3(point,down_x,down_z,up_y)
        if(a .lt. 0) then
          a=g_3(point,up_y,down_z,down_x)
          if(a .lt. 0) then
            a=g_3(point,up_y,down_x,down_z)
            if(a .lt. 0) then
              a=g_3(point,down_z,down_x,up_y)
              if(a .lt. 0) then
                a=g_3(point,down_z,up_y,down_x)
              end if
            end if
          end if
        end if
      end if
      adj(-1,1,-1)=a
c
      a=-1
      a=g_3(point,up_x,down_y,down_z)
      if(a .lt. 0) then
        a=g_3(point,up_x,down_z,down_y)
        if(a .lt. 0) then
          a=g_3(point,down_y,down_z,up_x)
          if(a .lt. 0) then
            a=g_3(point,down_y,up_x,down_z)
            if(a .lt. 0) then
              a=g_3(point,down_z,up_x,down_y)
              if(a .lt. 0) then
                a=g_3(point,down_z,down_y,up_x)
              end if
            end if
          end if
        end if
      end if
      adj(1,-1,-1)=a
c
      a=-1
      a=g_3(point,down_x,down_y,down_z)
      if(a .lt. 0) then
        a=g_3(point,down_x,down_z,down_y)
        if(a .lt. 0) then
          a=g_3(point,down_y,down_z,down_x)
          if(a .lt. 0) then
            a=g_3(point,down_y,down_x,down_z)
            if(a .lt. 0) then
              a=g_3(point,down_z,down_x,down_y)
              if(a .lt. 0) then
                a=g_3(point,down_z,down_y,down_x)
              end if
            end if
          end if
        end if
      end if
      adj(-1,-1,-1)=a
c
      return
      end
c
      integer function g_2(i,dir_1,dir_2)
c
      integer i,dir_1(*),dir_2(*)
c
      g_2=-1
c
      if(dir_1(i) .gt. 0) g_2=dir_2(dir_1(i))
c
      return
      end
c
      integer function g_3(i,dir_1,dir_2,dir_3)
c
      integer i,dir_1(*),dir_2(*),dir_3(*)
c
      g_3=-1
c
      if(dir_1(i) .gt. 0) then
         if(dir_2(dir_1(i)) .gt. 0) g_3=dir_3(dir_2(dir_1(i)))
      end if
c
      return
      end
