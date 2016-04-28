      subroutine find_uppers_hope(list_high,high_x,high_y,high_z,
     >  pos_high_x,pos_high_y,pos_high_z,
     >  high_point,hoc_daughter_points,next_daughter_points,
     >  pos_point_x,pos_point_y,pos_point_z,delta,
     >  point_up_x,point_up_y,point_up_z,point_pointer,points_maxx)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer points_maxx
      integer list_high,high_x(maxx),high_y(maxx),high_z(maxx)
      integer pos_high_x(maxx),pos_high_y(maxx),pos_high_z(maxx)
      integer high_point,hoc_daughter_points(maxx)
      integer next_daughter_points(points_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx),p_x,p_y,p_z,m,h_point
      integer p_p,point,p,point_pointer(points_maxx)
      integer n_x,n_y,n_z,pos,delta
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
c
      logical badd_x,badd_y,badd_z
c     
      point=hoc_daughter_points(high_point)
      do while(point .gt. 0)
        p_p=-point_pointer(point)
        p_x=-1
        p_y=-1
        p_z=-1
        m=(list_high-1)*27+p_p
c     
        h_point=high_x(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_x .lt. 0)
            if(pos_high_x(m) .eq. -point_pointer(p)) p_x=p              
            p=next_daughter_points(p)
          end do
          if(p_x .lt. 0) write(49,*)'x ',high_point,point,m,h_point,p_p,
     >      pos_high_x(m)
        end if
c     
        h_point=high_y(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_y .lt. 0)
            if(pos_high_y(m) .eq. -point_pointer(p)) p_y=p              
            p=next_daughter_points(p)
          end do
          if(p_y .lt. 0) write(49,*)'y ',high_point,point,m,h_point,p_p,
     >      pos_high_y(m)
        end if
c     
        h_point=high_z(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_z .lt. 0)
            if(pos_high_z(m) .eq. -point_pointer(p)) p_z=p       
            p=next_daughter_points(p)
          end do
          if(p_z .lt. 0) write(49,*)'z ',high_point,point,m,h_point,p_p,
     >      pos_high_z(m)
        end if
c     
        point_up_x(point)=p_x
        point_up_y(point)=p_y
        point_up_z(point)=p_z
c
        if(p_x .gt. 0) then
          badd_x=
     >      pos_point_x(p_x) .ne. pos_point_x(point)+delta .or.
     >      pos_point_y(p_x) .ne. pos_point_y(point) .or.
     >      pos_point_z(p_x) .ne. pos_point_z(point)
          if(badd_x) then
            write(49,*)'pos x ',pos_point_x(p_x),pos_point_x(point),
     >        pos_point_y(p_x),pos_point_y(point),
     >        pos_point_z(p_x),pos_point_z(point)
          end if
        end if
c
        if(p_y .gt. 0) then
          badd_y=
     >      pos_point_x(p_y) .ne. pos_point_x(point) .or.
     >      pos_point_y(p_y) .ne. pos_point_y(point)+delta .or.
     >      pos_point_z(p_y) .ne. pos_point_z(point)
          if(badd_y) then
            write(49,*)'pos y ',pos_point_x(p_y),pos_point_x(point),
     >        pos_point_y(p_y),pos_point_y(point),
     >        pos_point_z(p_y),pos_point_z(point)
          end if
        end if
c
        if(p_z .gt. 0) then
          badd_z=
     >      pos_point_x(p_z) .ne. pos_point_x(point) .or.
     >      pos_point_y(p_z) .ne. pos_point_y(point) .or.
     >      pos_point_z(p_z) .ne. pos_point_z(point)+delta
          if(badd_z) then
            write(49,*)'pos z ',pos_point_x(p_z),pos_point_x(point),
     >        pos_point_y(p_z),pos_point_y(point),
     >        pos_point_z(p_z),pos_point_z(point)
          end if
        end if
c     
        point=next_daughter_points(point)
      end do
      return
      end
