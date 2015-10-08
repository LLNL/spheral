      subroutine make_decisions_down(decisions_x,decisions_y,
     >  decisions_z,positions_x,positions_y,positions_z)
c     
      implicit none
c     
      logical decisions_x(27,27),decisions_y(27,27)
      logical decisions_z(27,27)
c     
      integer p_x_h,p_y_h,p_z_h,p_x_l,p_y_l,p_z_l
      integer pos_x_h,pos_y_h,pos_z_h,pos_x_l,pos_y_l,pos_z_l
      integer d_x,d_y,d_z,pos_h,pos_l,total_x,total_y,total_z
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
c
      data total_x,total_y,total_z/3*0/
c
      do p_z_h=1,3
       pos_z_h=(p_z_h-2)*2
       do p_y_h=1,3
        pos_y_h=(p_y_h-2)*2
        do p_x_h=1,3
         pos_x_h=(p_x_h-2)*2 
         pos_h=(p_z_h-1)*9+(p_y_h-1)*3+p_x_h
          do p_z_l=1,3
          pos_z_l=p_z_l-1
          d_z=pos_z_l-pos_z_h
          do p_y_l=1,3
           pos_y_l=p_y_l-1
           d_y=pos_y_l-pos_y_h
           do p_x_l=1,3
            pos_x_l=p_x_l-1
            pos_l=(p_z_l-1)*9+(p_y_l-1)*3+p_x_l
            d_x=pos_x_l-pos_x_h
c
            decisions_x(pos_l,pos_h)=
     >        d_x-1 .ge. 0 .and. d_x-1 .le. 2 .and.
     >        d_y .ge. 0 .and. d_y .le. 2 .and.
     >        d_z .ge. 0 .and. d_z .le. 2
            if(decisions_x(pos_l,pos_h))positions_x(pos_l,pos_h)=
     >        d_x+d_y*3+d_z*9
            if(decisions_x(pos_l,pos_h)) total_x=total_x+1
c
            decisions_y(pos_l,pos_h)=
     >        d_x .ge. 0 .and. d_x .le. 2 .and.
     >        d_y-1 .ge. 0 .and. d_y-1 .le. 2 .and.
     >        d_z .ge. 0 .and. d_z .le. 2
            if(decisions_y(pos_l,pos_h))positions_y(pos_l,pos_h)=
     >        (d_x+1)+(d_y-1)*3+d_z*9
            if(decisions_y(pos_l,pos_h)) total_y=total_y+1
c
            decisions_z(pos_l,pos_h)=
     >        d_x .ge. 0 .and. d_x .le. 2 .and.
     >        d_y .ge. 0 .and. d_y .le. 2 .and.
     >        d_z-1 .ge. 0 .and. d_z-1 .le. 2
            if(decisions_z(pos_l,pos_h))positions_z(pos_l,pos_h)=
     >        (d_x+1)+d_y*3+(d_z-1)*9
            if(decisions_z(pos_l,pos_h)) total_z=total_z+1
c
            print*,'decisions down ',pos_l,pos_h,
     >        positions_x(pos_l,pos_h),
     >        positions_y(pos_l,pos_h),positions_z(pos_l,pos_h),
     >        decisions_x(pos_l,pos_h),decisions_y(pos_l,pos_h),
     >        decisions_z(pos_l,pos_h)
           end do
          end do
         end do
        end do
       end do
      end do
c
      print*,'totals ',total_x,total_y,total_z
      return
      end
c
