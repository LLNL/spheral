      subroutine find_uppers(high_group,high_point,hoc_high_points,
     >  next_high_points,hoc_daughter_points,next_daughter_points,
     >  periodic,grid_multiply,d_point,point_up_x,point_up_y,point_up_z,
     >  pos_point_x,pos_point_y,pos_point_z,point_pointer,
     >  groups_maxx,points_maxx)
c     
      implicit none
c     
      integer groups_maxx,points_maxx
c     
      include 'maxx.inc'
c     
      integer p_1_h,p_2_h,p_3_h,point_0,p_x_0,p_y_0,p_z_0
      integer high_group,hoc_high_points(groups_maxx)
      integer next_high_points(points_maxx)
      integer hoc_daughter_points(points_maxx)
      integer next_daughter_points(points_maxx)
      integer d_point,pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer i,l,h_point,point
      integer point_pointer(points_maxx),d_point_2
      integer p_x,p_y,p_z,grid_multiply,dist_p
      integer high_point,p_1_0,p_2_0,p_3_0,p_1_x,p_2_y,p_3_z
      integer p_h(27),n_h,n_high,p_p,wrapping
      integer pointer_h(27),pointer_x,pointer_y,pointer_z
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
c     
      logical periodic,doit_x,doit_y,doit_z,first_time
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
c     
      dist_p(i,l)=min(abs(i),abs(i-l),abs(i+l))
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
        call make_decisions(decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
      end if
c     
      p_1_h=pos_point_x(high_point)
      p_2_h=pos_point_y(high_point)
      p_3_h=pos_point_z(high_point)
      d_point_2=d_point*2
      wrapping=grid_multiply/d_point_2
c*******
      if(periodic) then
        h_point=hoc_high_points(high_group)
        n_high=0
        do while(h_point .gt. 0)
          if(dist_p(pos_point_x(h_point)-p_1_h,grid_multiply) 
     >      .le. d_point_2) then
            if(dist_p(pos_point_y(h_point)-p_2_h,grid_multiply) 
     >        .le. d_point_2) then
              if(dist_p(pos_point_z(h_point)-p_3_h,grid_multiply) 
     >          .le. d_point_2) then
                n_high=n_high+1
                p_h(n_high)=h_point
                pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2+1
                pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2+1
                pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2+1
                call wrapper(pointer_x,wrapping)
                call wrapper(pointer_y,wrapping)
                call wrapper(pointer_z,wrapping)
                pointer_h(n_high)=
     >            1+pointer_x+pointer_y*3+pointer_z*9
              end if
            end if
          end if
          h_point=next_high_points(h_point)
        end do
      else
        h_point=hoc_high_points(high_group)
        n_high=0
        do while(h_point .gt. 0)
          if(abs(pos_point_x(h_point)-p_1_h) .le. d_point_2) then
            if(abs(pos_point_y(h_point)-p_2_h) .le. d_point_2) then
              if(abs(pos_point_z(h_point)-p_3_h) .le. d_point_2) 
     >          then
                n_high=n_high+1
                p_h(n_high)=h_point
                pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2+1
                pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2+1
                pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2+1
                pointer_h(n_high)=
     >            1+pointer_x+pointer_y*3+pointer_z*9
              end if
            end if
          end if
          h_point=next_high_points(h_point)
        end do
      end if
c**** 
c     
      point_0=hoc_daughter_points(high_point)
      do while(point_0 .gt. 0)
        p_x=-1
        p_y=-1
        p_z=-1
c     
        p_1_0=pos_point_x(point_0)
        p_2_0=pos_point_y(point_0)
        p_3_0=pos_point_z(point_0)
        p_1_x=p_1_0+d_point
        p_2_y=p_2_0+d_point
        p_3_z=p_3_0+d_point
c     
        p_p=-point_pointer(point_0)
c     
        n_h=1
        do while(n_h .le. n_high)
          h_point=p_h(n_h)
c     
          doit_x=decisions_x(p_p,pointer_h(n_h))
          if(doit_x) then
            p_x_0=-positions_x(p_p,pointer_h(n_h))
            point=hoc_daughter_points(h_point)
            do while(point .gt. 0 .and. p_x .lt. 0)
              if(p_x_0 .eq. point_pointer(point)) p_x=point
              point=next_daughter_points(point)
            end do
            point_up_x(point_0)=p_x
          end if
c     
          doit_y=decisions_y(p_p,pointer_h(n_h))
          if(doit_y) then
            p_y_0=-positions_y(p_p,pointer_h(n_h))
            point=hoc_daughter_points(h_point)
            do while(point .gt. 0 .and. p_y .lt. 0)
              if(p_y_0 .eq. point_pointer(point)) p_y=point
              point=next_daughter_points(point)
            end do
            point_up_y(point_0)=p_y
          end if
c     
          doit_z=decisions_z(p_p,pointer_h(n_h))
          if(doit_z) then
            p_z_0=-positions_z(p_p,pointer_h(n_h))
            point=hoc_daughter_points(h_point)
            do while(point .gt. 0 .and. p_z .lt. 0)
              if(p_z_0 .eq. point_pointer(point)) p_z=point
              point=next_daughter_points(point)
            end do
            point_up_z(point_0)=p_z
          end if
c     
          if(p_x .gt. 0 .and. p_y .gt. 0 .and. p_z .gt. 0) n_h=n_high
          n_h=n_h+1
        end do
        point_0=next_daughter_points(point_0)
      end do
c     
      return
      end
