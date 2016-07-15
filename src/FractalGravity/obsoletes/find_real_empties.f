      subroutine find_real_empties(high_group,high_point,
     >  hoc_high_points,next_high_points,periodic,
     >  grid_multiply,d_point,pos_point_x,pos_point_y,pos_point_z,
     >  empty,groups_maxx,points_maxx)
c     
      implicit none
c     
      integer groups_maxx,points_maxx
c     
      include 'maxx.inc'
c     
      integer p_1_h,p_2_h,p_3_h,n_x,n_y,n_z
      integer high_group,hoc_high_points(groups_maxx)
      integer next_high_points(points_maxx)
      integer d_point,pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer i,l,h_point,d_point_2,grid_multiply,dist_p
      integer high_point,wrapping,n_empties
      integer pointer_x,pointer_y,pointer_z
c     
      logical periodic,empty(-1:1,-1:1,-1:1)
c     
      dist_p(i,l)=min(abs(i),abs(i-l),abs(i+l))
c
c***  call WRAPPER
c***  wraparound
c     
      do n_z=-1,1
         do n_y=-1,1
            do n_x=-1,1
               empty(n_x,n_y,n_z)=.true.
            end do
         end do
      end do
c     
      p_1_h=pos_point_x(high_point)
      p_2_h=pos_point_y(high_point)
      p_3_h=pos_point_z(high_point)
      d_point_2=d_point*2
      wrapping=grid_multiply/d_point_2
c     
      if(periodic) then
         n_empties=0
         h_point=hoc_high_points(high_group)
         do while(h_point .gt. 0 .and. n_empties .lt. 27)
            if(dist_p(pos_point_x(h_point)-p_1_h,grid_multiply) 
     >        .le. d_point_2) then
               if(dist_p(pos_point_y(h_point)-p_2_h,grid_multiply) 
     >           .le. d_point_2) then
                  if(dist_p(pos_point_z(h_point)-p_3_h,grid_multiply) 
     >              .le. d_point_2) then
                     pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2+1
                     pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2+1
                     pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2+1
                     call wrapper(pointer_x,wrapping)
                     call wrapper(pointer_y,wrapping)
                     call wrapper(pointer_z,wrapping)
                     empty(pointer_x-1,pointer_y-1,pointer_z-1)=.false.
                     n_empties=n_empties+1
                  end if
               end if
            end if
            h_point=next_high_points(h_point)
         end do
      else
         n_empties=0
         h_point=hoc_high_points(high_group)
         do while(h_point .gt. 0 .and. n_empties .lt. 27)
            if(abs(pos_point_x(h_point)-p_1_h) .le. d_point_2) then
               if(abs(pos_point_y(h_point)-p_2_h) .le. d_point_2) then
                  if(abs(pos_point_z(h_point)-p_3_h) .le. d_point_2) 
     >              then
                     pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2
                     pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2
                     pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2
                     empty(pointer_x,pointer_y,pointer_z)=.false.
                     n_empties=n_empties+1
                  end if
               end if
            end if
            h_point=next_high_points(h_point)
         end do
      end if
c     
      return
      end
