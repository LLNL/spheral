      subroutine my_list_sort_3(hoc,next,x,y,z,up_x,point_pointer,
     >  delta,periodic,length,remove_duplicates,re_order,points_maxx)
c     
      implicit none
c     
      integer space,maxx_20
      parameter (space=10000)
      include 'maxx.inc'
      parameter (maxx_20=maxx/20)
c     
      logical periodic,remove_duplicates,re_order
      integer points_maxx
      integer hoc,next(points_maxx),x(points_maxx),y(points_maxx)
      integer z(points_maxx),up_x(points_maxx),delta,length
      integer hoc_z(space),next_z(maxx)
      integer hoc_y(space),next_y(maxx_20)
      integer hoc_x(space),next_x(maxx_20)
      integer up_s(maxx_20),pointer(maxx_20)
      integer x_s(maxx_20),y_s(maxx_20)
      integer n_x,n_y,n_z,n_tot_x,n_tot_y,n_tot_z,point
      integer point_pointer(points_maxx),hoc_a,next_a(maxx_20)
      integer hoc_space_tmp,next_space_tmp,hoc_s,p
      integer next_s(maxx_20)
c     
      hoc_space_tmp=space
      next_space_tmp=points_maxx
c     
      point=hoc
c     do while(point .gt. 0)
c     write(93,93)point,x(point),y(point),z(point)
c     93     format(' a ',10i8)
c     point=next(point)
c     end do
c     
      call my_list_sort_1(hoc,next,z,delta,periodic,
     >  hoc_z,next_z,n_tot_z,length,hoc_space_tmp,next_space_tmp)
c     
      if(re_order)hoc=-1
c     
      do n_z=n_tot_z,1,-1
c
        hoc_s=-1
        p=0
        point=hoc_z(n_z)
        do while(point .gt. 0)
          p=p+1
          x_s(p)=x(point)
          y_s(p)=y(point)
          pointer(p)=point
          next_s(p)=hoc_s
          hoc_s=p
          point=next_z(point)
        end do
c     
      p=hoc_s
      do while(p .gt. 0)
         point=pointer(p)
      write(93,94)n_z,p,point,x(point),x_s(p),y(point),y_s(p),z(point)
 94   format(' b ',10i8)
      p=next_s(p)
      end do
c     
        call my_list_sort_1(hoc_s,next_s,y_s,delta,periodic,
     >    hoc_y,next_y,n_tot_y,length,hoc_space_tmp,next_space_tmp)
c     
        do n_y=n_tot_y,1,-1
c     
           p=hoc_y(n_y)
           do while(p .gt. 0)
              write(93,95)n_z,n_y,p,x_s(p),y_s(p),z(pointer(p))
 95           format(' c ',10i8)
              p=next_y(p)
           end do
          call my_list_sort_1(hoc_y(n_y),next_y,x_s,delta,periodic,
     >      hoc_x,next_x,n_tot_x,length,hoc_space_tmp,next_space_tmp)
c     
          write(93,*)'n_tot_x= ',n_tot_x
          do n_x=n_tot_x,-1
             p=hoc_x(n_x)
             do while(p .gt. 0)
                write(93,96)n_z,n_y,n_x,p,x_s(p),y_s(p),z(pointer(p))
 96             format(' d ',10i8)
                p=next_x(p)
             end do
          end do
c
          if(remove_duplicates) call remove_dupes(hoc_x,next_x,
     >      n_tot_x,point_pointer)
c     
          hoc_a=-1
          do n_x=1,n_tot_x
            point=hoc_x(n_x)
            do while(point .gt. 0)
              next_a(point)=hoc_a
              hoc_a=point
              point=next_x(point)
            end do
          end do
c     
          call uppers(hoc_a,next_a,x_s,up_s,delta,length,periodic)
c     
          if(re_order) then
             p=hoc_a
             do while(p .gt. 0)
               point=pointer(p)
               up_x(point)=up_s(p)
               next(point)=hoc
               hoc=point
c                write(92,*)point,x(point),y(point),z(point),up_x(point),
c     >               point_pointer(point)
               p=next_a(p)
             end do
          end if
       end do
      end do
c     
      return
      end
