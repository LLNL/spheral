      subroutine my_list_sort_3(hoc,next,x,y,z,up_x,point_pointer,
     >  delta,periodic,length,remove_duplicates,re_order,points_maxx)
c     
      implicit none
c     
      integer space,maxx_20
      parameter (space=10000)
      include 'maxx.inc'
c     
      parameter (maxx_20=maxx/20)
c     
      logical periodic,remove_duplicates,re_order
      integer points_maxx
      integer hoc,next(points_maxx),x(points_maxx),y(points_maxx)
      integer z(points_maxx),up_x(points_maxx),delta,length
      integer hoc_z(space),next_z(maxx)
      integer hoc_y(space),next_y(maxx_20)
      integer hoc_x(space),next_x(maxx_20)
      integer x_a(maxx_20),y_a(maxx_20),up_a(maxx_20)
      integer p,pointer(-1:maxx_20),point_pointer_a(maxx_20)
      integer n_x,n_y,n_z,n_tot_x,n_tot_y,n_tot_z,point
      integer point_pointer(points_maxx),hoc_a,next_a(maxx_20)
      integer hoc_space_tmp,next_space_tmp
c     
      pointer(-1)=-1
      pointer(0)=0
      hoc_space_tmp=space
      next_space_tmp=maxx
c     
      call my_list_sort_1(hoc,next,z,delta,periodic,
     >  hoc_z,next_z,n_tot_z,length,hoc_space_tmp,next_space_tmp)
c     
      if(re_order)hoc=-1
c     
      do n_z=n_tot_z,1,-1
c     
        hoc_a=-1
        point=hoc_z(n_z)
        p=0
        do while(point .gt. 0)
          p=p+1
          next_a(p)=hoc_a
          hoc_a=p            
          x_a(p)=x(point)
          y_a(p)=y(point)
          pointer(p)=point
          point_pointer_a(p)=point_pointer(point)
c
          point=next_z(point)
        end do
c     
        next_space_tmp=maxx_20
c     
        call my_list_sort_1(hoc_a,next_a,y_a,delta,periodic,
     >    hoc_y,next_y,n_tot_y,length,hoc_space_tmp,next_space_tmp)
c     
        do n_y=n_tot_y,1,-1
c     
          next_space_tmp=maxx_20
c     
          call my_list_sort_1(hoc_y(n_y),next_y,x_a,delta,periodic,
     >      hoc_x,next_x,n_tot_x,length,hoc_space_tmp,next_space_tmp)
c     
          if(remove_duplicates) call remove_dupes(hoc_x,next_x,
     >      n_tot_x,point_pointer_a,hoc_space_tmp,next_space_tmp,
     >      points_maxx)
c     
          hoc_a=-1
          do n_x=1,n_tot_x
            point=hoc_x(n_x)
            do while(point .gt. 0)
              next_a(point)=hoc_a
              hoc_a=point
c
              point=next_x(point)
            end do
          end do
c     
          call uppers(hoc_a,next_a,x_a,up_a,delta,length,periodic,
     >      maxx_20,maxx_20)
c     
          p=hoc_a
          do while(p .gt. 0)
            point=pointer(p)
            up_x(point)=pointer(up_a(p))
            if(re_order) then
              next(point)=hoc
              hoc=point
            end if
            p=next_a(p)
          end do
        end do
      end do
c     
      return
      end
