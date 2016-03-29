      subroutine make_decisions_erika(decisions,positions)
c     
      implicit none
c     
      logical decisions(27,27)
c     
      integer p_x_h,p_y_h,p_z_h,p_x_l,p_y_l,p_z_l
      integer pos_x_h,pos_y_h,pos_z_h,pos_x_l,pos_y_l,pos_z_l
      integer d_x,d_y,d_z,pos_h,pos_l,total
      integer positions(27,27)
c     
      data total/0/
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
                        positions(pos_l,pos_h)=-1
                        decisions(pos_l,pos_h)=
     >                       d_x .ge. 0 .and. d_x .le. 2 .and.
     >                       d_y .ge. 0 .and. d_y .le. 2 .and.
     >                       d_z .ge. 0 .and. d_z .le. 2
c     
                        if(decisions(pos_l,pos_h))then
                           positions(pos_l,pos_h)=(d_x+1)+d_y*3+d_z*9
                           total=total+1
                        end if
c     
c                        print*,'decisions 0 ',pos_l,pos_h,
c     >                       positions(pos_l,pos_h),
c     >                       decisions(pos_l,pos_h)
                     end do
                  end do
               end do
            end do
         end do
      end do
c     
      print*,'totals 0 ',total
      return
      end
c
