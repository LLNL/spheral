      subroutine high_points(group,hoc_points,next_points,
     >     particles_at_point,minimum_number,really_high,n_highs,
     >     p_x,p_y,p_z,number_masks,level_mask,pos_mask,
     >     grid_length,level_max,level,periodic,shrink_mask,
     >     it_is_high,groups_maxx,points_maxx,high_points_maxx,
     >     memory_value)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx,high_points_maxx,n_highs
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer point,particles_at_point(points_maxx),minimum_number
      integer really_high(high_points_max_group),memory_value
c     
      integer p_x(points_maxx),p_y(points_maxx),p_z(points_maxx)
      integer level(groups_maxx),grid_multiply,number_masks
      integer level_mask(2000),grid_length,level_max
c     
      real pos_mask(10000),grid_multiplier,x,y,z,move
c     
      logical check_high,it_is_high(-1:points_maxx)
      logical check_mask,periodic,shrink_mask
c
      it_is_high(-1)=.false.
      it_is_high(0)=.false.
      if(shrink_mask) then
         move=(1.0-2.0**(-level(group)))/float(grid_length)
      else
         move=0.0
      end if
      grid_multiply=grid_length*2**level_max
      grid_multiplier=1.0/float(grid_multiply)
c     
      n_highs=0
      point=hoc_points(group)
      do while(point .gt. 0)
         it_is_high(point)=
     >        check_high(point,particles_at_point,minimum_number,
     >        points_maxx)
c     
         if(it_is_high(point)) then
            if(number_masks .gt. 0) then
               x=float(p_x(point))*grid_multiplier
               y=float(p_y(point))*grid_multiplier
               z=float(p_z(point))*grid_multiplier
               it_is_high(point)=check_mask(
     >              level(group),number_masks,
     >              level_mask,pos_mask,x,y,z,move,periodic)
c     
            end if
         end if
c     
         if(it_is_high(point)) then
            n_highs=n_highs+1
            if(n_highs .gt. high_points_maxx) then
               memory_value=5
               return
            end if
            really_high(n_highs)=point
         end if
         point=next_points(point)
      end do
c     
      return
      end
c     
      logical function check_mask(lev,number_masks,level_mask,
     >     mask,p_x,p_y,p_z,move,periodic)
c     
      implicit none
c     
      integer number_masks,i
      logical periodic
      integer lev,level_mask(2*number_masks),m
      real p_x,p_y,p_z,mask(6*number_masks)
      real x,move,dist_1,dx,dy,dz
c     
      dist_1(x)=min(abs(x),abs(x+1.0),abs(x-1.0))
c     
      check_mask=.false.
      m=-5
c     
      do i=1,number_masks
         m=m+6
         if(level_mask(2*i-1) .gt. lev) then
            if(periodic) then
               dx=dist_1(p_x-mask(m))
               dy=dist_1(p_y-mask(m+1))
               dz=dist_1(p_z-mask(m+2))
            else
               dx=abs(p_x-mask(m))
               dy=abs(p_y-mask(m+1))
               dz=abs(p_z-mask(m+2))
            end if
            if(level_mask(2*i) .eq. 0) then
               if(dx .le. mask(m+3)-move) then
                  if(dy .le. mask(m+4)-move) then
                     if(dz .le. mask(m+5)-move) then
                        check_mask=.true.
                        return
                     end if
                  end if
               end if
            else if(level_mask(2*i) .eq. 1) then
               check_mask=
     >              (dx/(mask(m+3)-move))**2+
     >              (dy/(mask(m+4)-move))**2+
     >              (dz/(mask(m+5)-move))**2 .lt. 1.0
               if(check_mask) return
            else
               stop ' die in masks'
            end if
         end if
      end do
c     
      return
      end
c
