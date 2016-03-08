      subroutine find_high_pairs(group,hoc_points,next_points,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  it_is_high,number_high_points,seq_to_list,high_pairs,
     >  list_high_1,list_high_2,debug,
     >  groups_maxx,points_maxx,
     >  points_maxx_group,particles_maxx,tweaks)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx,points_maxx_group,particles_maxx
      integer point,group,hoc_points(groups_maxx)
      integer next_points(points_maxx),tweaks
      integer n_h,number_high_points,high_pairs
      integer list_to_seq(maxx),seq_to_list(high_points_max_group)
      integer p_up_x,p_up_y,p_up_z,point_up_x(points_maxx)
      integer point_up_y(points_maxx),point_up_z(points_maxx)
      integer point_down_x(points_maxx)
      integer point_down_y(points_maxx),point_down_z(points_maxx)
      integer list_high_1(3*high_points_max_group)
      integer list_high_2(3*high_points_max_group)
      integer high_point,n
c     
      logical debug,it_is_high(-1:points_maxx)
c     
c***  a point is a high density point if 
c***  (particles_at_point(point) .ge. minimum_number)
c***  The quantity "minimum_number" is a user defined input parameter
c***  two high density points p_i and p_j are equivalent if
c***  p_up_x(p_i)=p_j or p_up_y(p_i)=p_j or p_up_z(p_i)=p_j.
c***  A group is an equivalence class of equivalent points
c***  We generate two arrays list_high_1(n) and list_high_2(n) so that they are
c***  the n'th pair of equivalent points. n=list_to_seq(point) and
c***  point=seq_to_list(n) point back and forth between the link list and the 
c***  sequential arrays
c     
c     
      if(debug) print*,'enter find high_pairs ',group
      point=hoc_points(group)
      number_high_points=0
      high_pairs=0
c     
      do while(point .gt. 0)
         if(it_is_high(point)) then
            number_high_points=number_high_points+1
c            print*,'too many ',point,number_high_points
            list_to_seq(point)=number_high_points
            seq_to_list(number_high_points)=point
            if(debug) then
               write(13,13)number_high_points,point,list_to_seq(point),
     >              seq_to_list(number_high_points)
            end if
         end if
         point=next_points(point)
      end do
c     
      if(number_high_points .eq. 0) return
c
      do n_h=1,number_high_points
         high_point=seq_to_list(n_h)
         p_up_x=point_up_x(high_point)
         p_up_y=point_up_y(high_point)
         p_up_z=point_up_z(high_point)
c     
         if(it_is_high(p_up_x))then
            high_pairs=high_pairs+1
            list_high_1(high_pairs)=n_h
            list_high_2(high_pairs)=list_to_seq(p_up_x)
         end if
c     
         if(it_is_high(p_up_y))then
            high_pairs=high_pairs+1
            list_high_1(high_pairs)=n_h
            list_high_2(high_pairs)=list_to_seq(p_up_y)
         end if
c     
         if(it_is_high(p_up_z))then
            high_pairs=high_pairs+1
            list_high_1(high_pairs)=n_h
            list_high_2(high_pairs)=list_to_seq(p_up_z)
         end if
c     
      end do
      if(debug)print*,'made it through'
      if(debug) then
         do n=1,high_pairs
            write(13,13)n,list_high_1(n),list_high_2(n),
     >           seq_to_list(list_high_1(n)),
     >           seq_to_list(list_high_2(n))
            if(list_high_1(n) .eq. list_high_2(n)) write(13,*)'one '
 13         format(7i8)
         end do
      end if
c     
      if(debug)print*,'exit find high_pairs'
      return
      end
c
