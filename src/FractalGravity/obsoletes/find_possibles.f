      subroutine find_possibles(high_group,new_group,hoc_high_points,
     >  next_high_points,pos,grid_length,level_max,level,
     >  hoc_possibles,next_possibles,groups_maxx,points_maxx)
c     
      implicit none
c
      include 'maxx.inc'
c
      integer groups_maxx,points_maxx
      integer number,grid_length,level(groups_maxx),high_group
      integer dist,n,hoc_possibles(maxx),h_point,new_group
      integer high_level
      integer hoc_high_points(points_maxx),next_possibles(maxx)
      integer pos(points_maxx),level_max,next_high_points(points_maxx)
c
      high_level=level(new_group)-1
      number=grid_length*2**high_level
      dist=2**(level_max-high_level)
c
      write(84,*)'possibles ',high_group,high_level,number
      do n=1,number
        hoc_possibles(n)=-1
      end do
c
      h_point=hoc_high_points(high_group)
      do while(h_point .gt. 0)
        n=pos(h_point)/dist+1
        next_possibles(h_point)=hoc_possibles(n)
        hoc_possibles(n)=h_point
        h_point=next_high_points(h_point)
      end do
c
      do n=1,number
        h_point=hoc_possibles(n)
        write(84,*)'well ',n,h_point
        do while(h_point .gt. 0)
c          write(84,*)n,h_point,pos(h_point)
          h_point=next_possibles(h_point)
        end do
      end do
c
      return
      end
