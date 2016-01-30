      subroutine find_go_ahead(go_ahead,ins,adj,decisions,positions,
     >  new_group)
c     
      implicit none
c     
      logical go_ahead(27),ins(27),pos,neg,decisions(27,27)
c     
      integer new_group
      integer adj(27)
      integer positions(27,27)
      integer p_l,p_h,i,hoc(27),next(27,27)
      logical first_time
      save hoc,next
c     
      pos(i)=i .gt. 0
      neg(i)=i .lt. 0
c     
      data first_time/.true./
c     
      if(first_time) then
        first_time=.false.
        do p_l=1,27
          hoc(p_l)=-1
          do p_h=1,27
            if(p_h .ne. 14) then
              if(decisions(p_l,p_h)) then
                if(positions(p_l,p_h) .lt. p_l) then
                  next(p_l,p_h)=hoc(p_l)
                  hoc(p_l)=p_h
                  print*,' decisions f ',p_l,p_h
                end if
              end if
            end if
          end do
        end do
      end if
c     
      do p_l=1,27
        ins(p_l)=.false.
        go_ahead(p_l)=.true.
        p_h=hoc(p_l)
        do while(p_h .gt. 0 .and. go_ahead(p_l))
          go_ahead(p_l)=neg(adj(p_h))
          p_h=next(p_l,p_h)
        end do
      end do
c     
      ins(1)=
     >  pos(adj(1)) .and. pos(adj(2)) .and. pos(adj(4)) 
     >  .and. pos(adj(5))
     >  .and. pos(adj(10)) .and. pos(adj(11)) .and. pos(adj(13))
      ins(2)=pos(adj(2)) .and. pos(adj(5)) .and. pos(adj(11))
      ins(4)=pos(adj(4)) .and. pos(adj(5)) .and. pos(adj(13))
      ins(5)=pos(adj(5))
      ins(10)=pos(adj(10)) .and. pos(adj(11)) .and. pos(adj(13))
      ins(11)=pos(adj(11))
      ins(13)=pos(adj(13))
      ins(14)=.true.
c     
      return
      end
