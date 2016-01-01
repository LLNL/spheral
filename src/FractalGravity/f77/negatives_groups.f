      subroutine negatives_groups(groups_maxx,
     >  hoc_points,mother_of_high_group,number_high_groups,
     >  next_groups,level,level_high,particles_in_group,
     >  points_in_group,hoc_high_groups,
     >  next_high_groups,mother_group,tweaks)
c     
      implicit none
c     
      integer groups_maxx,group,tweaks
      integer hoc_points(groups_maxx),mother_of_high_group(groups_maxx)
      integer number_high_groups(groups_maxx)
      integer next_groups(groups_maxx),level(groups_maxx)
      integer level_high(groups_maxx),particles_in_group(groups_maxx)
      integer hoc_high_groups(groups_maxx),next_high_groups(groups_maxx)
      integer mother_group(groups_maxx),points_in_group(groups_maxx)
c     
      do group=1,groups_maxx
       particles_in_group(group)=-group
       points_in_group(group)=-group
       hoc_points(group)=-group
       mother_of_high_group(group)=-group
       number_high_groups(group)=-group
       next_groups(group)=-group
       level(group)=-group
       level_high(group)=-group
       hoc_high_groups(group)=-group
       next_high_groups(group)=-group
       mother_group(group)=-group
      end do
c     
      return
      end
