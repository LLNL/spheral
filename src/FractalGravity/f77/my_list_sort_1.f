      subroutine my_list_sort_1(hoc,next,x,delta,periodic,hoc_1,next_1,
     >  n_tot_x,length,hoc_space_tmp,next_space_tmp)
c     
      implicit none
c     
      logical periodic,wrapped,first
c     
      integer hoc_space_tmp,next_space_tmp
      integer hoc,next(next_space_tmp),x(next_space_tmp)
      integer hoc_1(hoc_space_tmp),next_1(next_space_tmp)
      integer x_min,x_max,n_tot_x
      integer length,point,delta,n,x_min_a,x_min_b,x_max_a,x_max_b,half
c     
      point=hoc
      if(hoc .lt. 0) then
         n_tot_x=1
         hoc_1(1)=-1
         return
      end if
c
      x_min=2*length
      x_max=-1
      first=.false.
c
      do while(point .gt. 0)
         first=first .or. x(point) .eq. delta
         x_min=min(x_min,x(point))
         x_max=max(x_max,x(point))
         point=next(point)
      end  do
      wrapped=periodic .and. x_min .eq. 0 .and. first .and. 
     >     x_max .eq. length-delta
c     
      if(wrapped) then
        x_min_a=2*length
        x_min_b=2*length
        x_max_a=-1
        x_max_b=-1
        half=length/2
        point=hoc
        do while(point .gt. 0)
          if(x(point) .ge. half) then
            x_min_a=min(x_min_a,x(point))
            x_max_a=max(x_max_a,x(point))
          else
            x_min_b=min(x_min_b,x(point))
            x_max_b=max(x_max_b,x(point))
          end if
          point=next(point)
        end  do
        x_max=x_max_b
        x_min=x_min_a
        n_tot_x=(x_max-x_min+length)/delta+1
        do n=1,n_tot_x
          hoc_1(n)=-1
        end do
c     
        point=hoc
        do while(point .gt. 0)
          n=mod(x(point)-x_min+length,length)/delta+1
          if(n .le. 0 .or. n .gt. n_tot_x) then
            print*,point,x_min,x_max,length,n,n_tot_x
            stop 'baad in 1 '
          end if
          next_1(point)=hoc_1(n)
          hoc_1(n)=point
          point=next(point)
        end  do        
      else
c     
        n_tot_x=(x_max-x_min)/delta+1
        do n=1,n_tot_x
          hoc_1(n)=-1
        end do
c     
        point=hoc
        do while(point .gt. 0)
          n=(x(point)-x_min)/delta+1
          if(n .le. 0 .or. n .gt. n_tot_x) then
            print*,point,x_min,x_max,length,n,n_tot_x
            stop 'bad in 1 '
          end if
          next_1(point)=hoc_1(n)
          hoc_1(n)=point
          point=next(point)
        end  do
      end if
c     
      return
      end
