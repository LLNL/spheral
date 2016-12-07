      subroutine remove_duplicates(high_group,high_point,
     >  hoc_high_points,next_high_points,hoc_points,next_points,
     >  hoc_daughter_points,next_daughter_points,new_group,
     >  periodic,grid_multiply,d_point,hoc_backwards,next_backwards,
     >  pos_point_x,pos_point_y,pos_point_z,point_pointer,
     >  hoc_possibles,next_possibles,
     >  inside,groups_maxx,points_maxx)
c     
      implicit none
c     
      integer groups_maxx,points_maxx
c     
      include 'maxx.inc'
c     
      integer p_1_h,p_2_h,p_3_h,point_0,p_0
      integer high_group,hoc_high_points(groups_maxx)
      integer next_high_points(points_maxx)
      integer hoc_daughter_points(points_maxx)
      integer next_daughter_points(points_maxx)
      integer d_point,pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer i,l,h_point,point
      integer point_pointer(points_maxx),d_point_2
      integer grid_multiply,dist_p
      integer high_point,hoc_possibles(maxx),next_possibles(maxx)
      integer p_h(27),n_h,n_high,p_p,wrapping
      integer pointer_h(27),pointer_x,pointer_y,pointer_z
      integer positions_0(27,27),hoc_backwards
      integer next_backwards(maxx),up,down,point_tmp,n,hoc(3)
      integer hoc_points(groups_maxx),next_points(points_maxx),new_group
c     
      logical periodic
      logical doit_0,first_time
      logical decisions_0(27,27),inside(points_maxx)
c     
      dist_p(i,l)=min(abs(i),abs(i-l),abs(i+l))
c     
      data first_time/.true./
      data decisions_0/729*.false./
      data positions_0/729* -1000/
c     
      if(first_time) then
        first_time=.false.
        call make_decisions_0(decisions_0,positions_0)
      end if
c     
      p_1_h=pos_point_x(high_point)
      p_2_h=pos_point_y(high_point)
      p_3_h=pos_point_z(high_point)
      d_point_2=d_point*2
      wrapping=grid_multiply/d_point_2
c     
      n=p_1_h/d_point_2+1
      hoc(2)=hoc_possibles(n)
c     
      if(periodic) then
        if(n .eq. 1) then
          hoc(1)=hoc_possibles(wrapping)
        else
          hoc(1)=hoc_possibles(n-1)
        end if
        if(n .eq. wrapping) then
          hoc(3)=hoc_possibles(1)
        else
          hoc(3)=hoc_possibles(n+1)
        end if
      else
        hoc(1)=hoc_possibles(n-1)
        hoc(3)=hoc_possibles(n+1)
      end if
c     
c     
      if(periodic) then
        n_high=0
        do n=1,3
          h_point=hoc(n)
          do while(h_point .gt. 0)
            if(dist_p(pos_point_x(h_point)-p_1_h,grid_multiply) 
     >        .le. d_point_2) then
              if(dist_p(pos_point_y(h_point)-p_2_h,grid_multiply) 
     >          .le. d_point_2) then
                if(dist_p(pos_point_z(h_point)-p_3_h,grid_multiply) 
     >            .le. d_point_2) then
                  n_high=n_high+1
                  p_h(n_high)=h_point
                  pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2+1
                  pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2+1
                  pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2+1
                  call wrapper(pointer_x,wrapping)
                  call wrapper(pointer_y,wrapping)
                  call wrapper(pointer_z,wrapping)
                  pointer_h(n_high)=1+pointer_x+pointer_y*3+pointer_z*9
                end if
              end if
            end if
            h_point=next_possibles(h_point)
          end do
        end do
      else
        n_high=0
        do n=1,3
          h_point=hoc(n)
          do while(h_point .gt. 0)
            if(abs(pos_point_x(h_point)-p_1_h) .le. d_point_2) then
              if(abs(pos_point_y(h_point)-p_2_h) .le. d_point_2) then
                if(abs(pos_point_z(h_point)-p_3_h) .le. d_point_2) then
                  n_high=n_high+1
                  p_h(n_high)=h_point
                  pointer_x=(pos_point_x(h_point)-p_1_h)/d_point_2+1
                  pointer_y=(pos_point_y(h_point)-p_2_h)/d_point_2+1
                  pointer_z=(pos_point_z(h_point)-p_3_h)/d_point_2+1
                  pointer_h(n_high)=1+pointer_x+pointer_y*3+pointer_z*9
                end if
              end if
            end if
            h_point=next_possibles(h_point)
          end do
        end do
      end if
c     
      point_0=hoc_daughter_points(high_point)
      do while(point_0 .gt. 0)
        if(.not. inside(point_0)) then
c     
          p_p=-point_pointer(point_0)
c     
          n_h=1
          do while(n_h .le. n_high)
            h_point=p_h(n_h)
c     
            doit_0=decisions_0(p_p,pointer_h(n_h)) 
     >        .and. high_point .ne. h_point
            if(doit_0) then
              point=hoc_daughter_points(h_point)
              point_tmp=-1
              p_0=positions_0(p_p,pointer_h(n_h))
              do while(point .gt. 0)
                if(p_0 .eq. -point_pointer(point)) then
                  if(-point_pointer(point) .gt. p_p) then
                    write(81,*)'removed ',point_0,p_p,
     >                point,-point_pointer(point)
                    if(point_tmp .gt. 0) then
                      next_daughter_points(point_tmp)=
     >                  next_daughter_points(point)
                    else
                      hoc_daughter_points(h_point)=
     >                  next_daughter_points(
     >                  hoc_daughter_points(h_point))
                    end if
c     
                    up=next_points(point)
                    down=next_backwards(point)
                    if(down .gt. 0) then
                      next_points(down)=up
                    else
                      hoc_points(new_group)=up
                    end if
                    if(up .gt. 0) then
                      next_backwards(up)=down
                    else
                      hoc_backwards=down
                    end if
                  end if
                end if
                point_tmp=point
                point=next_daughter_points(point)
              end do
            end if
c     
            n_h=n_h+1
          end do
        end if
        point_0=next_daughter_points(point_0)
      end do
c     
      return
      end
