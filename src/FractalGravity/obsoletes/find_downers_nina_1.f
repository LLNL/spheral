      subroutine find_downers_nina(adj,periodic,grid_multiply,
     >  high_x,high_y,high_z,
     >  pos_high_x,pos_high_y,pos_high_z,
     >  high_point,hoc_daughter_points,next_daughter_points,
     >  pos_point_x,pos_point_y,pos_point_z,delta,
     >  point_down_x,point_down_y,point_down_z,point_pointer,
     >  hoc_high,next_high_points,debug,high_points_maxx,points_maxx)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer high_points_maxx,points_maxx
      integer high_x(27),high_y(27),high_z(27),adj(27)
      integer pos_high_x(27),pos_high_y(27),pos_high_z(27)
      integer high_point,hoc_daughter_points(points_maxx)
      integer next_daughter_points(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx),p_x,p_y,p_z,m
      integer h_point,grid_multiply
      integer p_p,point,p,point_pointer(points_maxx)
      integer delta,hoc_high,next_high_points(maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),i,j,k
c     
      logical badd_x,badd_y,badd_z,first_time,must_exist(3,27)
      logical periodic,problem_x,problem_y,problem_z,dist_p,debug
c     
      dist_p(i,j,k)=min(abs(i-j),abs(i-j-k),abs(i-j+k)) .ne. 0
c     
      save must_exist
      data first_time/.true./
c     
      if(first_time) then
        first_time=.false.
        do p_z=1,3
          do p_y=1,3
            do p_x=1,3
              p=(p_z-1)*9+(p_y-1)*3+p_x
              must_exist(1,p)=p_x .gt. 1
              must_exist(2,p)=p_y .gt. 1
              must_exist(3,p)=p_z .gt. 1
            end do
          end do
        end do
      end if
c     
      point=hoc_daughter_points(high_point)
      do while(point .gt. 0)
        p_p=-point_pointer(point)
        p_x=-1
        p_y=-1
        p_z=-1
        m=p_p
c     
        h_point=high_x(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_x .lt. 0)
            if(pos_high_x(m) .eq. -point_pointer(p)) p_x=p              
            p=next_daughter_points(p)
          end do
          if(p_x .lt. 0) then
            call emergency_search_f(periodic,grid_multiply,
     >        hoc_high,next_high_points,
     >        hoc_daughter_points,next_daughter_points,pos_point_x,
     >        pos_point_y,pos_point_z,point,-delta,0,0,delta,p_x,
     >        high_points_maxx,points_maxx)
            write(49,*)'xd ',high_point,point,h_point,p_p,
     >        pos_high_x(m),p_x
          end if
        end if
c     
        h_point=high_y(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_y .lt. 0)
            if(pos_high_y(m) .eq. -point_pointer(p)) p_y=p              
            p=next_daughter_points(p)
          end do
          if(p_y .lt. 0) then
            call emergency_search_f(periodic,grid_multiply,
     >        hoc_high,next_high_points,
     >        hoc_daughter_points,next_daughter_points,pos_point_x,
     >        pos_point_y,pos_point_z,point,0,-delta,0,delta,p_y,
     >        high_points_maxx,points_maxx)
            write(49,*)'yd ',high_point,point,h_point,p_p,
     >        pos_high_y(m),p_y
          end if
        end if
c     
        h_point=high_z(m)
        if(h_point .gt. 0) then
          p=hoc_daughter_points(h_point)
          do while(p .gt. 0 .and. p_z .lt. 0)
            if(pos_high_z(m) .eq. -point_pointer(p)) p_z=p       
            p=next_daughter_points(p)
          end do
          if(p_z .lt. 0) then
            call emergency_search_f(periodic,grid_multiply,
     >        hoc_high,next_high_points,
     >        hoc_daughter_points,next_daughter_points,pos_point_x,
     >        pos_point_y,pos_point_z,point,0,0,-delta,delta,p_z,
     >        high_points_maxx,points_maxx)
            write(49,*)'zd ',high_point,point,h_point,p_p,
     >        pos_high_z(m),p_z
          end if
        end if
c     
        point_down_x(point)=p_x
        point_down_y(point)=p_y
        point_down_z(point)=p_z
c     
        if(debug) then
          badd_x=.false.
          badd_y=.false.
          badd_z=.false.
c     
          if(periodic) then
            if(p_x .gt. 0) then
              badd_x=
     >          dist_p(pos_point_x(p_x),pos_point_x(point)-delta,
     >          grid_multiply) .or.
     >          dist_p(pos_point_y(p_x),pos_point_y(point),
     >          grid_multiply) .or.
     >          dist_p(pos_point_z(p_x),pos_point_z(point),
     >          grid_multiply)
            end if
c     
            if(p_y .gt. 0) then
              badd_y=
     >          dist_p(pos_point_x(p_y),pos_point_x(point),
     >          grid_multiply) .or.
     >          dist_p(pos_point_y(p_y),pos_point_y(point)-delta,
     >          grid_multiply) .or.
     >          dist_p(pos_point_z(p_y),pos_point_z(point),
     >          grid_multiply)
            end if
c     
            if(p_z .gt. 0) then
              badd_z=
     >          dist_p(pos_point_x(p_z),pos_point_x(point),
     >          grid_multiply) .or.
     >          dist_p(pos_point_y(p_z),pos_point_y(point),
     >          grid_multiply) .or.
     >          dist_p(pos_point_z(p_z),pos_point_z(point)-delta,
     >          grid_multiply)
            end if
c     
          else
            if(p_x .gt. 0) then
              badd_x=
     >          pos_point_x(p_x) .ne. pos_point_x(point)-delta .or.
     >          pos_point_y(p_x) .ne. pos_point_y(point) .or.
     >          pos_point_z(p_x) .ne. pos_point_z(point)
            end if
c     
            if(p_y .gt. 0) then
              badd_y=
     >          pos_point_x(p_y) .ne. pos_point_x(point) .or.
     >          pos_point_y(p_y) .ne. pos_point_y(point)-delta .or.
     >          pos_point_z(p_y) .ne. pos_point_z(point)
            end if
c     
            if(p_z .gt. 0) then
              badd_z=
     >          pos_point_x(p_z) .ne. pos_point_x(point) .or.
     >          pos_point_y(p_z) .ne. pos_point_y(point) .or.
     >          pos_point_z(p_z) .ne. pos_point_z(point)-delta
            end if
c     
          end if
c     
          if(badd_x) then
            write(49,*)'pos x d ',badd_x,delta,
     >        pos_point_x(p_x),pos_point_x(point),
     >        pos_point_y(p_x),pos_point_y(point),
     >        pos_point_z(p_x),pos_point_z(point)
          end if
c     
          if(badd_y) then
            write(49,*)'pos y d ',badd_y,delta,
     >        pos_point_x(p_y),pos_point_x(point),
     >        pos_point_y(p_y),pos_point_y(point),
     >        pos_point_z(p_y),pos_point_z(point)
          end if
c     
          if(badd_z) then
            write(49,*)'pos z d ',badd_z,delta,
     >        pos_point_x(p_z),pos_point_x(point),
     >        pos_point_y(p_z),pos_point_y(point),
     >        pos_point_z(p_z),pos_point_z(point)
          end if
c     
          problem_x=must_exist(1,m) .and. point_down_x(point) .le. 0
          problem_y=must_exist(2,m) .and. point_down_y(point) .le. 0
          problem_z=must_exist(3,m) .and. point_down_z(point) .le. 0
c     
          if(problem_x .or. problem_y .or. problem_z) then
            write(59,*)high_point,m,point,point_down_x(point),
     >        point_down_y(point),point_down_z(point),
     >        problem_x,problem_y,problem_z
          end if
        end if
        point=next_daughter_points(point)
      end do
      return
      end
