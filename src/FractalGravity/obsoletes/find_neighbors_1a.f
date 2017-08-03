      subroutine find_neighbors(point,adj,
     >  up_x,up_y,up_z,down_x,down_y,down_z,delta,
     >  edge,spacing,box_length,hoc_highs,next_highs,
     >  pos_point_x,pos_point_y,pos_point_z,
     >  points_maxx)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer points_maxx
      integer point,i,j,k,up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx)
      integer down_y(points_maxx),down_z(points_maxx)
      integer adj(-1:1,-1:1,-1:1),g_2,g_3
      integer n,n_x,n_y,n_z,p_x,p_y,p_z,h_point,delta
c     
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer edge(6),spacing(3),box_length(3),hoc_highs(length_maxx)
      integer next_highs(length_maxx)
c     
      logical dont_bother(-1:1,-1:1,-1:1)
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
      adj(-1,0,0)=down_x(point)
      adj(0,-1,0)=down_y(point)
      adj(0,0,-1)=down_z(point)
c     
      adj(-1,-1,0)=g_2(point,down_x,down_y)
      adj(1,-1,0)=g_2(point,up_x,down_y)
      adj(-1,1,0)=g_2(point,down_x,up_y)
      adj(1,1,0)=g_2(point,up_x,up_y)
      adj(1,0,1)=g_2(point,up_x,up_z)
      adj(0,1,1)=g_2(point,up_y,up_z)
      adj(-1,0,1)=g_2(point,down_x,up_z)
      adj(0,-1,1)=g_2(point,down_y,up_z)
      adj(1,0,-1)=g_2(point,up_x,down_z)
      adj(0,1,-1)=g_2(point,up_y,down_z)
      adj(-1,0,-1)=g_2(point,down_x,down_z)
      adj(0,-1,-1)=g_2(point,down_y,down_z)
c     
      adj(1,1,1)=g_3(point,up_x,up_y,up_z)
      adj(-1,1,1)=g_3(point,down_x,up_y,up_z)
      adj(1,-1,1)=g_3(point,up_x,down_y,up_z)
      adj(-1,-1,1)=g_3(point,down_x,down_y,up_z)
      adj(1,1,-1)=g_3(point,up_x,up_y,down_z)
      adj(-1,1,-1)=g_3(point,down_x,up_y,down_z)
      adj(1,-1,-1)=g_3(point,up_x,down_y,down_z)
      adj(-1,-1,-1)=g_3(point,down_x,down_y,down_z)
c     
      do k=-1,1
        do j=-1,1
          do i=-1,1
            dont_bother(i,j,k)=adj(i,j,k) .gt. 0
          end do
        end do
      end do
c     
c     dont_bother(0,0,0)=.true.
      do i=-1,1,2
        dont_bother(i,0,0)=.true.
        dont_bother(0,i,0)=.true.
        dont_bother(0,0,i)=.true.
      end do
c     
      do j=-1,1,2
        do i=-1,1,2
          dont_bother(i,j,0)=dont_bother(i,j,0) .or.
     >      adj(i,0,0) .gt. 0 .or. adj(0,j,0) .gt. 0
          dont_bother(i,0,j)=dont_bother(i,0,j) .or.
     >      adj(i,0,0) .gt. 0 .or. adj(0,0,j) .gt. 0
          dont_bother(0,i,j)=dont_bother(0,i,j) .or.
     >      adj(0,i,0) .gt. 0 .or. adj(0,0,j) .gt. 0
        end do
      end do
c     
      do k=-1,1,2
        do j=-1,1,2
          do i=-1,1,2
            dont_bother(i,j,k)=dont_bother(i,j,k) .or.
     >        adj(0,j,k) .gt. 0 .or. adj(i,0,k) .gt. 0 .or.
     >        adj(i,j,0) .gt. 0
          end do
        end do
      end do
c     
      do k=-1,1
        do j=-1,1
          do i=-1,1
            if(.not. dont_bother(i,j,k)) then
              p_x=pos_point_x(point)+delta*i
              p_y=pos_point_y(point)+delta*j
              p_z=pos_point_z(point)+delta*k
              n_x=(p_x-edge(1))/spacing(1)
              n_y=(p_y-edge(3))/spacing(2)
              n_z=(p_z-edge(5))/spacing(3)
              n=1+n_x+box_length(1)*(n_y+box_length(2)*n_z)
              h_point=hoc_highs(n)
              do while(h_point .gt. 0 .and. adj(i,j,k) .lt. 0)
                if(p_z .eq. pos_point_z(h_point)) then
                  if(p_y .eq. pos_point_y(h_point)) then
                    if(p_x .eq. pos_point_x(h_point)) then
                      adj(i,j,k)=h_point
                    end if
                  end if
                end if
                h_point=next_highs(h_point)
              end do
            end if
          end do
        end do
      end do      
c     
      return
      end
c
      integer function g_2(i,dir_1,dir_2)
c
      integer i,dir_1(*),dir_2(*)
c
      g_2=-1
      if(dir_1(i) .gt. 0) g_2=dir_2(dir_1(i))
      if(g_2 .lt. 0) then
        if(dir_2(i) .gt. 0) g_2=dir_1(dir_2(i))
      end if
c
      return
      end
c
      integer function g_3(i,dir_1,dir_2,dir_3)
c
      integer i,dir_1(*),dir_2(*),dir_3(*),g_2
c
      g_3=g_2(i,dir_1,dir_2)
      if(g_3 .gt. 0) g_3=dir_3(g_3)
c
      if(g_3 .lt. 0) then
        g_3=g_2(i,dir_2,dir_3)
        if(g_3 .gt. 0) g_3=dir_1(g_3)
        if(g_3 .lt. 0) then
          g_3=g_2(i,dir_1,dir_3)
          if(g_3 .gt. 0) g_3=dir_2(g_3)
        end if
      end if
c
      return
      end
