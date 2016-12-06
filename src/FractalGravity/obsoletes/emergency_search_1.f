      subroutine emergency_search(periodic,grid_multiply,
     >  hoc_high,next_high_points,
     >  hoc_daughter_points,next_daughter_points,pos_point_x,
     >  pos_point_y,pos_point_z,point_0,delta_1,delta_2,delta_3,
     >  delta,point_u,high_points_maxx,points_maxx)
c     
      implicit none
c     
      integer high_points_maxx,points_maxx
      integer hoc_high,next_high_points(high_points_maxx)
      integer hoc_daughter_points(high_points_maxx)
      integer next_daughter_points(points_maxx),i,j,k,m,delta
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),point_0,delta_1,delta_2,delta_3
      integer p_x,p_y,p_z,point,point_u,high_point,grid_multiply
      integer p_x_d,p_y_d,p_z_d
c     
      logical periodic,dist_p,dist_p_d
c     
      dist_p(i,j,k)=min(abs(i-j),abs(i-j+k),abs(i-j-k)) .eq. 0
c     
      dist_p_d(i,j,k,m)=min(abs(i-j),abs(i-j+k),abs(i-j-k)) .le. m
c     
      p_x=pos_point_x(point_0)+delta_1
      p_y=pos_point_y(point_0)+delta_2
      p_z=pos_point_z(point_0)+delta_3
      p_x_d=p_x-delta
      p_y_d=p_y-delta
      p_z_d=p_z-delta
      point_u=-1
c     
      if(periodic) then
        high_point=hoc_high
        do while(high_point .gt. 0 .and. point_u .lt. 0)
          if(dist_p_d(p_x_d,pos_point_x(high_point),
     >      grid_multiply,delta)) then
            if(dist_p_d(p_y_d,pos_point_y(high_point),
     >        grid_multiply,delta)) then
              if(dist_p_d(p_z_d,pos_point_z(high_point),
     >          grid_multiply,delta)) then
                point=hoc_daughter_points(high_point)
                do while(point .gt. 0 .and. point_u .lt. 0)
                  if(dist_p(p_x,pos_point_x(point),grid_multiply)) 
     >              then
                    if(dist_p(p_y,pos_point_y(point),grid_multiply)) 
     >                then
                      if(dist_p(p_z,pos_point_z(point),grid_multiply)) 
     >                  then
                        point_u=point
                      end if
                    end if
                  end if
                  point=next_daughter_points(point)
                end do
              end if
            end if
          end if
          high_point=next_high_points(high_point)
        end do
      else
        high_point=hoc_high
        do while(high_point .gt. 0 .and. point_u .lt. 0)
          if(abs(p_x_d-pos_point_x(high_point)) .le. delta) then
            if(abs(p_y_d-pos_point_z(high_point)) .le. delta) then
              if(abs(p_z_d-pos_point_z(high_point)) .le. delta) then
                point=hoc_daughter_points(high_point)
                do while(point .gt. 0 .and. point_u .lt. 0)
                  if(p_x .eq. pos_point_x(point)) then
                    if(p_y .eq. pos_point_y(point)) then
                      if(p_z .eq. pos_point_z(point)) then
                        point_u=point
                      end if
                    end if
                  end if
                  point=next_daughter_points(point)
                end do
              end if
            end if
          end if
          high_point=next_high_points(high_point)
        end do
      end if
c     
      return
      end
