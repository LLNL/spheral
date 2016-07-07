      subroutine make_decisions_nina(decisions_xm,decisions_ym,
     >  decisions_zm,positions_xm,positions_ym,positions_zm)
c     
      implicit none
c     
      logical decisions_xm(27,27),decisions_ym(27,27)
      logical decisions_zm(27,27)
c     
      integer p_x_h,p_y_h,p_z_h,p_x_l,p_y_l,p_z_l
      integer pos_x_h,pos_y_h,pos_z_h,pos_x_l,pos_y_l,pos_z_l
      integer d_x,d_y,d_z,pos_h,pos_l,total_x,total_y,total_z
      integer positions_xm(27,27),positions_ym(27,27)
      integer positions_zm(27,27)
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
                  decisions_xm(pos_l,pos_h)=
     >              d_x-1 .ge. 0 .and. d_x-1 .le. 2 .and.
     >              d_y .ge. 0 .and. d_y .le. 2 .and.
     >              d_z .ge. 0 .and. d_z .le. 2
                  positions_xm(pos_l,pos_h)=-1
                  if(decisions_xm(pos_l,pos_h))
     >              positions_xm(pos_l,pos_h)=d_x+d_y*3+d_z*9
                  if(decisions_xm(pos_l,pos_h)) total_x=total_x+1
c     
                  positions_ym(pos_l,pos_h)=-1
                  decisions_ym(pos_l,pos_h)=
     >              d_x .ge. 0 .and. d_x .le. 2 .and.
     >              d_y-1 .ge. 0 .and. d_y-1 .le. 2 .and.
     >              d_z .ge. 0 .and. d_z .le. 2
                  if(decisions_ym(pos_l,pos_h))
     >              positions_ym(pos_l,pos_h)=(d_x+1)+(d_y-1)*3+d_z*9
                  if(decisions_ym(pos_l,pos_h)) total_y=total_y+1
c     
                  positions_zm(pos_l,pos_h)=-1
                  decisions_zm(pos_l,pos_h)=
     >              d_x .ge. 0 .and. d_x .le. 2 .and.
     >              d_y .ge. 0 .and. d_y .le. 2 .and.
     >              d_z-1 .ge. 0 .and. d_z-1 .le. 2
                  if(decisions_zm(pos_l,pos_h))
     >              positions_zm(pos_l,pos_h)=(d_x+1)+d_y*3+(d_z-1)*9
                  if(decisions_zm(pos_l,pos_h)) total_z=total_z+1
c     
                  print*,'decisions down ',pos_l,pos_h,
     >              positions_xm(pos_l,pos_h),
     >              positions_ym(pos_l,pos_h),positions_zm(pos_l,pos_h),
     >              decisions_xm(pos_l,pos_h),decisions_ym(pos_l,pos_h),
     >              decisions_zm(pos_l,pos_h)
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
