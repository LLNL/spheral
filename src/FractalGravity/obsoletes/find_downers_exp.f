      subroutine find_downers(high_point,hoc_daughter_points,
     >  point_down_x,point_down_y,point_down_z,
     >  next_daughter_points,point_pointer,adj,points_maxx)
c     
      implicit none
c     
      integer points_maxx
c     
      integer point_0,p_x_0,p_y_0,p_z_0
      integer hoc_daughter_points(points_maxx)
      integer next_daughter_points(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx)
      integer h_point,point
      integer point_pointer(points_maxx)
      integer p_x,p_y,p_z
      integer high_point,p_p,h_pointer
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
      integer adj(-1:1,-1:1,-1:1),n_x,n_y,n_z
c     
      logical doit_x,doit_y,doit_z,first_time
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
c     
      data first_time/.true./
      data decisions_x/729*.false./
      data decisions_y/729*.false./
      data decisions_z/729*.false./
      data positions_x/729* -1000/
      data positions_y/729* -1000/
      data positions_z/729* -1000/
c     
      if(first_time) then
        first_time=.false.
        call make_decisions_down(decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
      end if
c     
      point_0=hoc_daughter_points(high_point)
      do while(point_0 .gt. 0)
        p_x=-1
        p_y=-1
        p_z=-1
c     
        p_p=-point_pointer(point_0)
c     
        do n_z=-1,1
          do n_y=-1,1
            do n_x=-1,1
c     
              h_point=adj(n_x,n_y,n_z)
              if(h_point .gt. 0) then
c     
                h_pointer=n_z*9+n_y*3+n_x+14
                doit_x=decisions_x(p_p,h_pointer)
                if(doit_x) then
                  p_x_0=-positions_x(p_p,h_pointer)
                  point=hoc_daughter_points(h_point)
                  do while(point .gt. 0 .and. p_x .lt. 0)
                    if(p_x_0 .eq. point_pointer(point)) p_x=point
                    point=next_daughter_points(point)
                  end do
                  point_down_x(point_0)=p_x
                end if
c     
                doit_y=decisions_y(p_p,h_pointer)
                if(doit_y) then
                  p_y_0=-positions_y(p_p,h_pointer)
                  point=hoc_daughter_points(h_point)
                  do while(point .gt. 0 .and. p_y .lt. 0)
                    if(p_y_0 .eq. point_pointer(point)) p_y=point
                    point=next_daughter_points(point)
                  end do
                  point_down_y(point_0)=p_y
                end if
c     
                doit_z=decisions_z(p_p,h_pointer)
                if(doit_z) then
                  p_z_0=-positions_z(p_p,h_pointer)
                  point=hoc_daughter_points(h_point)
                  do while(point .gt. 0 .and. p_z .lt. 0)
                    if(p_z_0 .eq. point_pointer(point)) p_z=point
                    point=next_daughter_points(point)
                  end do
                  point_down_z(point_0)=p_z
                end if
              end if
            end do
          end do
        end do
        point_0=next_daughter_points(point_0)
      end do
c     
      return
      end
c        if(p_p .eq. 1) then
c          p_x=point_0+1
c          p_y=point_0+2
c          p_z=point_0+4
c        else if(p_p .eq. 2) then
c          p_y=point_0+2
c          p_z=point_0+4
c        else if(p_p .eq. 3) then
c          p_y=point_0+1
c          p_z=point_0+2
c        else if(p_p .eq. 4) then
c          p_x=point_0+1
c          p_z=point_0+4
c        else if(p_p .eq. 5) then
c          p_z=point_0+4
c        else if(p_p .eq. 6) then
c          p_z=point_0+2
c        else if(p_p .eq. 7) then
c          p_x=point_0+1
c          p_z=point_0+2
c        else if(p_p .eq. 8) then
c          p_z=point_0+2
c        else if(p_p .eq. 9) then
c          p_z=point_0+1
c        else if(p_p .eq. 10) then
c          p_x=point_0+1
c          p_y=point_0+2
c        else if(p_p .eq. 11) then
c          p_y=point_0+2
c        else if(p_p .eq. 12) then
c          p_y=point_0+1
c        else if(p_p .eq. 13) then
c          p_x=point_0+1
c        else if(p_p .eq. 16) then
c          p_x=point_0+1
c        else if(p_p .eq. 19) then
c          p_x=point_0+1
c          p_y=point_0+2
c        else if(p_p .eq. 20) then
c          p_y=point_0+2
c        else if(p_p .eq. 21) then
c          p_y=point_0+1
c        else if(p_p .eq. 22) then
c          p_x=point_0+1
c        else if(p_p .eq. 25) then
c          p_x=point_0+1
c        end if
