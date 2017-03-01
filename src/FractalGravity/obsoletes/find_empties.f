      subroutine find_empties(point,empty,
     >  up_x,up_y,up_z,down_x,down_y,down_z)
c     
      implicit none
c     
      integer point,i,j,k,up_x(*),up_y(*),up_z(*)
      integer down_x(*),down_y(*),down_z(*)
c     
      logical empty(-1:1,-1:1,-1:1),f_1,f_2,f_3
c     
      f_1(i)=i .lt. 0
c     
      do k=-1,1
       do j=-1,1
        do i=-1,1
         empty(i,j,k)=.true.
        end do
       end do
      end do
c     
      empty(0,0,0)=.false.
c     
      empty(1,0,0)=f_1(up_x(point))
      empty(0,1,0)=f_1(up_y(point))
      empty(0,0,1)=f_1(up_z(point))
c     
      empty(-1,0,0)=f_1(down_x(point))
      empty(0,-1,0)=f_1(down_y(point))
      empty(0,0,-1)=f_1(down_z(point))
c     
      empty(-1,-1,0)=
     >  f_2(point,down_x,down_y) .and.
     >  f_2(point,down_y,down_x)
      empty(1,-1,0)=
     >  f_2(point,up_x,down_y) .and.
     >  f_2(point,down_y,up_x)
      empty(-1,1,0)=
     >  f_2(point,down_x,up_y) .and.
     >  f_2(point,up_y,down_x)
      empty(1,1,0)=
     >  f_2(point,up_x,up_y) .and.
     >  f_2(point,up_y,up_x)
c     
      empty(1,0,1)=
     >  f_2(point,up_x,up_z) .and.
     >  f_2(point,up_z,up_x)
      empty(0,1,1)=
     >  f_2(point,up_y,up_z) .and.
     >  f_2(point,up_z,up_y)
      empty(-1,0,1)=
     >  f_2(point,down_x,up_z) .and.
     >  f_2(point,up_z,down_x)
      empty(0,-1,1)=
     >  f_2(point,down_y,up_z) .and.
     >  f_2(point,up_z,down_y)
c     
      empty(1,0,-1)=
     >  f_2(point,up_x,down_z) .and.
     >  f_2(point,down_z,up_x)
      empty(0,1,-1)=
     >  f_2(point,up_y,down_z) .and.
     >  f_2(point,down_z,up_y)
      empty(-1,0,-1)=
     >  f_2(point,down_x,down_z) .and.
     >  f_2(point,down_z,down_x)
      empty(0,-1,-1)=
     >  f_2(point,down_y,down_z) .and.
     >  f_2(point,down_z,down_y)
c     
      empty(1,1,1)=
     >  f_3(point,up_x,up_y,up_z) .and.
     >  f_3(point,up_x,up_z,up_y) .and.
     >  f_3(point,up_y,up_z,up_x) .and.
     >  f_3(point,up_y,up_x,up_z) .and.
     >  f_3(point,up_z,up_x,up_y) .and.
     >  f_3(point,up_z,up_y,up_x)
      empty(-1,1,1)=
     >  f_3(point,down_x,up_y,up_z) .and.
     >  f_3(point,down_x,up_z,up_y) .and.
     >  f_3(point,up_y,up_z,down_x) .and.
     >  f_3(point,up_y,down_x,up_z) .and.
     >  f_3(point,up_z,down_x,up_y) .and.
     >  f_3(point,up_z,up_y,down_x)
      empty(1,-1,1)=
     >  f_3(point,up_x,down_y,up_z) .and.
     >  f_3(point,up_x,up_z,down_y) .and.
     >  f_3(point,down_y,up_z,up_x) .and.
     >  f_3(point,down_y,up_x,up_z) .and.
     >  f_3(point,up_z,up_x,down_y) .and.
     >  f_3(point,up_z,down_y,up_x)
      empty(-1,-1,1)=
     >  f_3(point,down_x,down_y,up_z) .and.
     >  f_3(point,down_x,up_z,down_y) .and.
     >  f_3(point,down_y,up_z,down_x) .and.
     >  f_3(point,down_y,down_x,up_z) .and.
     >  f_3(point,up_z,down_x,down_y) .and.
     >  f_3(point,up_z,down_y,down_x)
c     
      empty(1,1,-1)=
     >  f_3(point,up_x,up_y,down_z) .and.
     >  f_3(point,up_x,down_z,up_y) .and.
     >  f_3(point,up_y,down_z,up_x) .and.
     >  f_3(point,up_y,up_x,down_z) .and.
     >  f_3(point,down_z,up_x,up_y) .and.
     >  f_3(point,down_z,up_y,up_x)
      empty(-1,1,-1)=
     >  f_3(point,down_x,up_y,down_z) .and.
     >  f_3(point,down_x,down_z,up_y) .and.
     >  f_3(point,up_y,down_z,down_x) .and.
     >  f_3(point,up_y,down_x,down_z) .and.
     >  f_3(point,down_z,down_x,up_y) .and.
     >  f_3(point,down_z,up_y,down_x)
      empty(1,-1,-1)=
     >  f_3(point,up_x,down_y,down_z) .and.
     >  f_3(point,up_x,down_z,down_y) .and.
     >  f_3(point,down_y,down_z,up_x) .and.
     >  f_3(point,down_y,up_x,down_z) .and.
     >  f_3(point,down_z,up_x,down_y) .and.
     >  f_3(point,down_z,down_y,up_x)
      empty(-1,-1,-1)=
     >  f_3(point,down_x,down_y,down_z) .and.
     >  f_3(point,down_x,down_z,down_y) .and.
     >  f_3(point,down_y,down_z,down_x) .and.
     >  f_3(point,down_y,down_x,down_z) .and.
     >  f_3(point,down_z,down_x,down_y) .and.
     >  f_3(point,down_z,down_y,down_x)
      return
      end
c
      logical function f_2(i,dir_1,dir_2)
c
      integer i,dir_1(*),dir_2(*)
c
      f_2=.true.
c
      if(dir_1(i) .gt. 0) f_2=dir_2(dir_1(i)) .lt. 0
c
      return
      end
c
      logical function f_3(i,dir_1,dir_2,dir_3)
c
      integer i,dir_1(*),dir_2(*),dir_3(*)
c
      f_3=.true.
c
      if(dir_1(i) .gt. 0) then
         if(dir_2(dir_1(i)) .gt. 0) f_3=dir_3(dir_2(dir_1(i))) .lt. 0
      end if
c
      return
      end
