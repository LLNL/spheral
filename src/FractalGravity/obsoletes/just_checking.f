      subroutine just_checking(adj,pos_point_x,pos_point_y,pos_point_z,
     >  hoc_daughter_points,next_daughter_points,delta,points_maxx)
c     
      implicit none
c     
      integer points_maxx
      integer adj(27),pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),hoc_daughter_points(points_maxx)
      integer next_daughter_points(points_maxx),delta,high_point
      integer p,n_x,n_y,n_z,p_h
      integer h_point,point,n_points(27),p_x,p_y,p_z
c     
      high_point=adj(14)
      p=0
      do n_z=0,2
        do n_y=0,2
          do n_x=0,2
            p=p+1
            n_points(p)=0
            p_x=pos_point_x(high_point)+delta*n_x
            p_y=pos_point_y(high_point)+delta*n_y
            p_z=pos_point_z(high_point)+delta*n_z
            do p_h=1,27
              h_point=adj(p_h)
              if(h_point .gt. 0) then
                point=hoc_daughter_points(h_point)
                do while(point .gt. 0)
                  if(pos_point_x(point) .eq. p_x) then
                    if(pos_point_y(point) .eq. p_y) then
                      if(pos_point_z(point) .eq. p_z) then
                        n_points(p)=n_points(p)+1
                      end if
                    end if
                  end if   
                  point=next_daughter_points(point)
                end do
              end if
            end do
c     
            if(n_points(p) .ne. 1) then
              write(58,*)'wrong ',high_point,p,n_points(p),p_x,p_y,p_z
c              problems=problems+1
c              p_pointer(problems)=p
c              what_problem(problems)=n_points(p)
            end if
          end do
        end do
      end do
c     
      return
      end
