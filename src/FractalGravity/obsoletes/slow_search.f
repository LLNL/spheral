      subroutine slow_search(adj,periodic,grid_multiply,
     >  hoc_daughter_points,next_daughter_points,pos_point_x,
     >  pos_point_y,pos_point_z,point_0,delta_1,delta_2,delta_3,
     >  point_u,high_points_maxx,points_maxx)
c     
      implicit none
c
      logical dist_p,periodic
c
      integer high_points_maxx,points_maxx
      integer adj(27),h_p,grid_multiply,i,j,k
      integer hoc_daughter_points(high_points_maxx)
      integer next_daughter_points(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),point_0,delta_1,delta_2,delta_3
      integer p_x,p_y,p_z,point,point_u,high_point
c
      dist_p(i,j,k)=min(abs(i-j),abs(i-j+k),abs(i-j-k)) .eq. 0
c
c      if(points_maxx .gt.0) return
c     
      p_x=pos_point_x(point_0)+delta_1
      p_y=pos_point_y(point_0)+delta_2
      p_z=pos_point_z(point_0)+delta_3
      point_u=-1
c     
      if(periodic) then
        do h_p=1,27
          high_point=adj(h_p)
          if (high_point .gt. 0 .and. point_u .lt. 0) then
            point=hoc_daughter_points(high_point)
            do while(point .gt. 0 .and. point_u .lt. 0)
              if(dist_p(p_x,pos_point_x(point),grid_multiply)) then
                if(dist_p(p_z,pos_point_y(point),grid_multiply)) then
                  if(dist_p(p_z,pos_point_z(point),grid_multiply)) then
                    point_u=point
                    write(49,*)'found it'
                  end if
                end if
              end if
              point=next_daughter_points(point)
            end do
          end if
        end do
      else
        do h_p=1,27
          high_point=adj(h_p)
          if (high_point .gt. 0 .and. point_u .lt. 0) then
            point=hoc_daughter_points(high_point)
            do while(point .gt. 0 .and. point_u .lt. 0)
              if(p_x .eq. pos_point_x(point)) then
                if(p_y .eq. pos_point_y(point)) then
                  if(p_z .eq. pos_point_z(point)) then
                    point_u=point
                    write(49,*)'found it'
                  end if
                end if
              end if
              point=next_daughter_points(point)
            end do
          end if
        end do
      end if
c     
      return
      end
