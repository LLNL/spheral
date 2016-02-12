      subroutine high_groups(highest_used_high_group,
     >  number_high_groups,mother_group_current,
     >  hoc_high_groups,next_high_groups,
     >  number_high_points,hoc_high_points,
     >  next_high_points,
     >  mother_of_high_group,seq_to_list,group_tmp,
     >  it_is_high,level,level_high,
     >  point_up_x,point_up_y,point_up_z,point_down_x,point_down_y,
     >  point_down_z,
     >  debug,groups_maxx,points_maxx,high_points_maxx,tweaks)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx,high_points_maxx,tweaks
      integer highest_used_high_group,number_high_groups(groups_maxx)
      integer mother_group_current,hoc_high_groups(groups_maxx)
      integer next_high_groups(groups_maxx)
      integer seq_to_list(high_points_max_group)
      integer family(high_points_max_group)
      integer number_high_points,hoc_high_points(groups_maxx)
      integer next_high_points(high_points_maxx)
      integer mother_of_high_group(groups_maxx),n
      integer high_group,group_tmp(high_points_max_group)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx),high_point
      integer level(groups_maxx),level_high(groups_maxx)
c     
      logical debug,it_is_high(-1:points_maxx)
c     
c     
c***  first assign each point to a group
c***  point(n) belongs to the same group as point(group_tmp(n))
c***  point=seq_to_list(n) points back to the single link list
c***  we generate individual link-list for each group 
c***  high_point=hoc_high_points(high_group) and 
c***  high_point=next_high_points(high_point)
c***  high_point_(up,down)_(x,y)(high_point)=point_(up,down)_(x,y)(high_point)
c***  if particles_at_point(point_(up,down)_(x,y)(high_point)) .
c***  ge. minimum_number
c***  else high_point_(up,down)_(x,y)(high_point)=-1
c     
      if(debug) print*,'enter high_groups'
c     
      high_group=highest_used_high_group
      number_high_groups(mother_group_current)=0
      hoc_high_groups(mother_group_current)=-1
c     
      do n=1,number_high_points
        if(group_tmp(n) .eq. n) then
          high_point=seq_to_list(n)
          number_high_groups(mother_group_current)=
     >      number_high_groups(mother_group_current)+1
          high_group=high_group+1
          level_high(high_group)=level(mother_group_current)
          hoc_high_points(high_group)=high_point
          next_high_points(high_point)=-1
          family(n)=high_group
          mother_of_high_group(high_group)=mother_group_current
          next_high_groups(high_group)=
     >      hoc_high_groups(mother_group_current)
          hoc_high_groups(mother_group_current)=high_group
        end if
      end do
      highest_used_high_group=high_group
c     
      do n=1,number_high_points
        if(group_tmp(n) .ne. n) then
          high_group=family(group_tmp(n))
          high_point=seq_to_list(n)
          next_high_points(high_point)=hoc_high_points(high_group)
          hoc_high_points(high_group)=high_point
        end if
      end do
c     
      if(debug) then
        high_group=hoc_high_groups(mother_group_current)
        do while(high_group .gt. 0)
          high_point=hoc_high_points(high_group)
          do while(high_point .gt. 0)
            if(debug)write(13,*)'high points ',high_group,high_point
            high_point=next_high_points(high_point)
          end do
          high_group=next_high_groups(high_group)
        end do       
      end if
c     
      if(debug) print*,'exit high_groups'
      if(debug) write(32,*)'exit high_groups'
      return
      end
