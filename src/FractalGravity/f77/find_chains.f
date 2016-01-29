      subroutine find_chains(hoc_groups,next_groups,mother_group,
     >  level,number_high_groups,number_chains,hoc_groups_chain,
     >  next_groups_chain,done_group,chain_started,chain_finished,
     >  debug,groups_maxx,tweaks)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,number_chains,group,chain,tweaks
      integer hoc_groups,next_groups(groups_maxx)
      integer number_high_groups(groups_maxx),level(groups_maxx)
      integer hoc_groups_chain(groups_maxx),mother_group(groups_maxx)
      integer next_groups_chain(groups_maxx)
c     
      logical done_group(groups_maxx),chain_started(groups_maxx)
      logical chain_finished(groups_maxx),used(g_maxx),debug
c     
c      write(41,*)'find chains started'
      number_chains=0
      group=hoc_groups
      do while(group .gt. 0)
        if(debug)write(41,*)group,number_chains
        if(group .gt. 1) then
          used(group)=.false.
          done_group(group)=.false.
          if(number_high_groups(group) .le. 0) then
            number_chains=number_chains+1
            hoc_groups_chain(number_chains)=group
            next_groups_chain(group)=-1
            used(group)=.true.
            chain_started(number_chains)=.false.
            chain_finished(number_chains)=.false.
          end if
        end if
        group=next_groups(group)
      end do
c     
      group=1
      done_group(group)=.true.
      used(group)=.true.
c     
      if(number_chains .le. 0) return
c     
      do chain=1,number_chains
        group=hoc_groups_chain(chain)
        do while(.not. used(mother_group(group)))
          used(mother_group(group))=.true.
          hoc_groups_chain(chain)=mother_group(group)
          next_groups_chain(mother_group(group))=group
          group=mother_group(group)
        end do
      end do
c     
      if(debug) then
        do chain=1,number_chains
          group=hoc_groups_chain(chain)
          do while(group .gt. 0) 
            write(37,*)'chain ',chain,group,level(group),
     >        mother_group(group)
            group=next_groups_chain(group)
          end do
        end do
      end if
      return
      end
