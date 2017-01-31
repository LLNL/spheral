      subroutine test_pos(high_group,hoc_high_points,next_high_points,
     >  new_group,hoc_points,next_points,pos_point_x,pos_point_y,
     >  pos_point_z,point_pointer,level_max,level,periodic,
     >  grid_multiply)
c     
      implicit none
c     
      logical found,periodic
      integer high_group,hoc_high_points(*),next_high_points(*)
      integer new_group,hoc_points(*),next_points(*)
      integer pos_point_x(*),pos_point_y(*),pos_point_z(*)
      integer point_pointer(*),level_max,level(*)
      integer delta,h_point,i,j,k,p_x,p_y,p_z,point,grid_multiply
c     
      delta=2**(level_max-level(new_group))
      h_point=hoc_high_points(high_group)
      print*,'grid_multiply= ',grid_multiply
      do while(h_point .gt. 0)
        do k=0,2
          do j=0,2
            do i=0,2
              p_x=pos_point_x(h_point)+i*delta
              p_y=pos_point_y(h_point)+j*delta
              p_z=pos_point_z(h_point)+k*delta
              if(periodic) then
                p_x=mod(p_x+grid_multiply,grid_multiply)
                p_y=mod(p_y+grid_multiply,grid_multiply)
                p_z=mod(p_z+grid_multiply,grid_multiply)
              end if
              point=hoc_points(new_group)
              found=.false.
              do while(point .gt. 0)
                if(pos_point_x(point) .eq. p_x) then
                  if(pos_point_y(point) .eq. p_y) then
                    if(pos_point_z(point) .eq. p_z) then
                      if(found) then
                        write(44,*)'ohoh ',new_group,h_point,
     >                    p_x,p_y,p_z,
     >                    i,j,k,point_pointer(point),point
                      else
                        found=.true.
                        write(44,*)'good ',new_group,h_point,
     >                    p_x,p_y,p_z,
     >                    i,j,k,point_pointer(point),point
                      end if
                    end if
                  end if
                end if
                point=next_points(point)
              end do
              if(.not. found) then
                write(44,*)'badd ',new_group,h_point,p_x,p_y,p_z,
     >            i,j,k
              end if
            end do
          end do
        end do
        h_point=next_high_points(h_point)
      end do
      return
      end
