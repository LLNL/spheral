      subroutine look_where(adj,list_high,
     >  high_x,high_y,high_z,
     >  pos_high_x,pos_high_y,pos_high_z,
     >  decisions_x,decisions_y,decisions_z,
     >  positions_x,positions_y,positions_z)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer adj(-1:1,-1:1,-1:1),a_1,a_2,a_3,a_x,a_y,a_z
      integer list_high,high_x(maxx),high_y(maxx),high_z(maxx)
      integer pos_high_x(maxx),pos_high_y(maxx),pos_high_z(maxx)
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
      integer pos_l,n_x_l,n_y_l,n_z_l,h_x,h_y,h_z,p_x,p_y,p_z
      integer p_x_min,p_y_min,p_z_min,n_x_h,n_y_h,n_z_h,pos_h,m
      integer tmp_x,tmp_y,tmp_z
c     
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
c
      a_1(m)=mod(m-1,3)-1
      a_2(m)=mod((m-1)/3,3)-1
      a_3(m)=(m-1)/9-1
c     
      do pos_l=1,27
        n_x_l=mod(pos_l-1,3)-1
        n_y_l=mod((pos_l-1)/3,3)-1
        n_z_l=(pos_l-1)/9-1
        h_x=-1
        h_y=-1
        h_z=-1
        p_x=-1
        p_y=-1
        p_z=-1
        tmp_x=-1
        tmp_y=-1
        tmp_z=-1
        p_x_min=28
        p_y_min=28
        p_z_min=28
        do pos_h=1,27
          n_x_h=mod(pos_h-1,3)-1
          n_y_h=mod((pos_h-1)/3,3)-1
          n_z_h=(pos_h-1)/9-1
          if(adj(n_x_h,n_y_h,n_z_h) .gt. 0) then
            if(decisions_x(pos_l,pos_h)) then
              if(positions_x(pos_l,pos_h) .lt. p_x_min) then
                h_x=pos_h
                p_x=positions_x(pos_l,pos_h)
                p_x_min=p_x
              end if
            end if
            if(decisions_y(pos_l,pos_h)) then
              if(positions_y(pos_l,pos_h) .lt. p_y_min) then
                h_y=pos_h
                p_y=positions_y(pos_l,pos_h)
                p_y_min=p_y
              end if
            end if
            if(decisions_z(pos_l,pos_h)) then
              if(positions_z(pos_l,pos_h) .lt. p_z_min) then
                h_z=pos_h
                p_z=positions_z(pos_l,pos_h)
                p_z_min=p_z
              end if
            end if
          end if
c     
        end do
        m=(list_high-1)*27+pos_l
        if(h_x .gt. 0) then
          tmp_x=h_x
          a_x=a_1(h_x)
          a_y=a_2(h_x)
          a_z=a_3(h_x)
          h_x=adj(a_x,a_y,a_z)
        end if
        high_x(m)=h_x
        pos_high_x(m)=p_x
        if(h_y .gt. 0) then
          tmp_y=h_y
          a_x=a_1(h_y)
          a_y=a_2(h_y)
          a_z=a_3(h_y)
          h_y=adj(a_x,a_y,a_z)
        end if
        high_y(m)=h_y
        pos_high_y(m)=p_y
        if(h_z .gt. 0) then
          tmp_z=h_z
          a_x=a_1(h_z)
          a_y=a_2(h_z)
          a_z=a_3(h_z)
          h_z=adj(a_x,a_y,a_z)
        end if
        high_z(m)=h_z
        pos_high_z(m)=p_z
        write(42,42)list_high,m,pos_l,h_x,p_x,h_y,p_y,h_z,p_z,
     >    tmp_x,tmp_y,tmp_z
 42     format('list ',12i8)
      end do

      return
      end
