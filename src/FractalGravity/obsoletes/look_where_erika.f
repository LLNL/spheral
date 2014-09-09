      subroutine look_where_erika(adj,high_point,
     >  high_x,high_y,high_z,
     >  pos_high_x,pos_high_y,pos_high_z,
     >  decisions_x,decisions_y,decisions_z,
     >  positions_x,positions_y,positions_z)
c     
      implicit none
c     
      integer adj(27)
      integer high_x(27),high_y(27),high_z(27)
      integer pos_high_x(27),pos_high_y(27),pos_high_z(27)
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
      integer pos_l,h_x,h_y,h_z,p_x,p_y,p_z
      integer p_x_min,p_y_min,p_z_min,pos_h
      integer tmp_x,tmp_y,tmp_z,high_point
c     
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
c     
      do pos_l=1,27
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
c     
        do pos_h=1,27
          if(adj(pos_h) .gt. 0) then
            if(decisions_x(pos_l,pos_h)) then
              if(positions_x(pos_l,pos_h) .lt. p_x_min) then
                h_x=pos_h
                p_x=positions_x(pos_l,pos_h)
                p_x_min=p_x
              end if
            end if
c     
            if(decisions_y(pos_l,pos_h)) then
              if(positions_y(pos_l,pos_h) .lt. p_y_min) then
                h_y=pos_h
                p_y=positions_y(pos_l,pos_h)
                p_y_min=p_y
              end if
            end if
c     
            if(decisions_z(pos_l,pos_h)) then
              if(positions_z(pos_l,pos_h) .lt. p_z_min) then
                h_z=pos_h
                p_z=positions_z(pos_l,pos_h)
                p_z_min=p_z
              end if
            end if
          end if
        end do
c     
        if(h_x .gt. 0) then
          tmp_x=h_x
          h_x=adj(h_x)
        end if
        high_x(pos_l)=h_x
        pos_high_x(pos_l)=p_x
c     
        if(h_y .gt. 0) then
          tmp_y=h_y
          h_y=adj(h_y)
        end if
        high_y(pos_l)=h_y
        pos_high_y(pos_l)=p_y
c     
        if(h_z .gt. 0) then
          tmp_z=h_z
          h_z=adj(h_z)
        end if
        high_z(pos_l)=h_z
        pos_high_z(pos_l)=p_z
c     
c        write(42,42)high_point,pos_l,
c     >    h_x,p_x,h_y,p_y,h_z,p_z,
c     >    tmp_x,tmp_y,tmp_z
c 42     format('list ',14i8)
      end do

      return
      end
