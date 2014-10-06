      subroutine remove_dupes(hoc,next,n_tot,pointer,
     >  hoc_space_tmp,next_space_tmp,points_maxx)
c
      implicit none
c
      integer hoc_space_tmp,next_space_tmp,points_maxx,n_tot
      integer hoc(hoc_space_tmp),next(next_space_tmp)
      integer pointer(points_maxx)
      integer n,point,point_0,p_0
c
      do n=1,n_tot
         point=hoc(n)
         if(point .gt. 0) then
            point_0=point
            p_0=pointer(point)
            point=next(point)
            do while(point .gt. 0)
               if(pointer(point) .gt. p_0) then
                  point_0=point
                  p_0=pointer(point)
               end if
               point=next(point)
            end do
            hoc(n)=point_0
            next(point_0)=-1
         end if
      end do
      return
      end
