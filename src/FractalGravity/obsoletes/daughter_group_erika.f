      subroutine daughter_group(new_group,hoc_points,
     >  next_points,hoc_high_points,next_high_points,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  hoc_particles,next_particles,pos_x,pos_y,
     >  pos_z,mother_group,high_group,inside,periodic,
     >  level,pos_point_x,pos_point_y,pos_point_z,
     >  point_pointer,level_max,high_point_up_x,
     >  high_point_up_y,high_point_up_z,
     >  particles_at_point,particles_in_group,
     >  points_in_group,high_point_down_x,
     >  high_point_down_y,high_point_down_z,
     >  particle_pointer,
     >  lowest_point_in_group,highest_point_in_group,
     >  highest_used_point,highest_used_particle,
     >  grid_length,debug,groups_maxx,points_maxx,
     >  points_maxx_group,high_points_maxx,
     >  particles_maxx,particles_real,memory_value,tweaks)
c     
      implicit none
c
c     
      integer groups_maxx,points_maxx,points_maxx_group,high_points_maxx
      integer particles_maxx,particles_real,tweaks
c     
      logical go_ahead(27),ins(27)
      logical periodic,inside(points_maxx),debug
      logical first_time
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
      logical decisions_xm(27,27),decisions_ym(27,27)
      logical decisions_zm(27,27)
      logical decisions(27,27)
c     
      include 'maxx.inc'
c     
      integer mother_group(groups_maxx),new_group
      integer d_point,level_max,level(groups_maxx),grid_length
      integer highest_used_particle,point,highest_used_point
      integer hoc_points(groups_maxx),high_point,delta
      integer hoc_high_points(groups_maxx),high_group,hoc_high
      integer point_down_x(points_maxx)
      integer point_down_y(points_maxx),point_down_z(points_maxx)
      integer n_p
      integer next_points(points_maxx),pos_point_x(points_maxx)
      integer pos_point_y(points_maxx),pos_point_z(points_maxx)
      integer hoc_particles(points_maxx),particles_at_point(points_maxx)
      integer point_pointer(points_maxx),particles_in_group(groups_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx),particle,new_particle
      integer next_particles(particles_maxx)
      integer next_daughter_points(maxx)
      integer points_in_group(groups_maxx)
      integer high_point_up_x(high_points_maxx)
      integer high_point_down_x(high_points_maxx)
      integer high_point_up_y(high_points_maxx)
      integer high_point_down_y(high_points_maxx)
      integer high_point_up_z(high_points_maxx)
      integer high_point_down_z(high_points_maxx)
      integer hoc_daughter_points(maxx),grid_multiply
      integer particle_pointer(particles_maxx)
      integer next_high_points(high_points_maxx)
      integer lowest_point_in_group(groups_maxx)
      integer highest_point_in_group(groups_maxx)
      integer memory_value
      integer adj(-1:1,-1:1,-1:1),total_points
      integer positions(27,27),p,p_x,p_y,p_z,point_tmp(27)
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
      integer positions_xm(27,27),positions_ym(27,27)
      integer positions_zm(27,27)
c     
      integer high_x(maxx),high_y(maxx),high_z(maxx)
      integer pos_high_x(maxx),pos_high_y(maxx),pos_high_z(maxx)
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),d_point_inv
c     
      real t_1,t_2,t_3,t_4,t_5
c     
      data first_time/.true./
      data decisions/729*.false./
      data decisions_x/729*.false./
      data decisions_y/729*.false./
      data decisions_z/729*.false./
      data decisions_xm/729*.false./
      data decisions_ym/729*.false./
      data decisions_zm/729*.false./
      data positions/729*-1000/
      data positions_x/729*-1000/
      data positions_y/729*-1000/
      data positions_z/729*-1000/
      data positions_xm/729*-1000/
      data positions_ym/729*-1000/
      data positions_zm/729*-1000/
c     
c     
      call cpu_time(t_1)
c     
c     print*,' I am here in daughter group'
c     
      if(first_time) then
        first_time=.false.
        call make_decisions_erika(decisions,positions)
        call make_decisions(decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
        call make_decisions_nina(decisions_xm,decisions_ym,decisions_zm,
     >    positions_xm,positions_ym,positions_zm)
      end if
c     
      memory_value=0
      level(new_group)=level(mother_group(new_group))+1
      d_point=2**(level_max-level(new_group))
      d_point_inv=1.0/float(d_point)
      grid_multiply=grid_length*2**level_max
      new_particle=highest_used_particle
      point=highest_used_point
      hoc_points(new_group)=-1
      points_in_group(new_group)=0
      particles_in_group(new_group)=0
c     
      high_point=hoc_high_points(high_group)
      hoc_high=high_point
c     
      do while(high_point .gt. 0)
c     
        if(point+27 .ge. points_maxx) then
          memory_value=1
          print*,'memory value= ',memory_value
          print*,'exceeded points in daughter ',point+27
          return
        end if
        if(point+27 .ge. maxx) then
          memory_value=3
          print*,'memory value= ',memory_value
          print*,'exceeded next_daughter in daughter ',point+27
          return
        end if
c
        call find_neighbors(high_point,adj,
     >    high_point_up_x,high_point_up_y,high_point_up_z,
     >    high_point_down_x,high_point_down_y,high_point_down_z,
     >    points_maxx)
c
        call find_go_ahead(go_ahead,ins,adj,decisions,positions,
     >    new_group)
c     
        hoc_daughter_points(high_point)=-1
c
        do p=1,27
          if(go_ahead(p)) then
            p_x=mod((p-1),3)
            p_y=mod((p-1)/3,3)
            p_z=(p-1)/9
c     
            point=point+1
            next_points(point)=hoc_points(new_group)
            hoc_points(new_group)=point
            pos_point_x(point)=pos_point_x(high_point)+
     >        p_x*d_point
            pos_point_y(point)=pos_point_y(high_point)+
     >        p_y*d_point
            pos_point_z(point)=pos_point_z(high_point)+
     >        p_z*d_point
            hoc_particles(point)=-1
            particles_at_point(point)=0
            next_daughter_points(point)=
     >        hoc_daughter_points(high_point)
            hoc_daughter_points(high_point)=point
            point_pointer(point)=-p
            point_tmp(p)=point
            points_in_group(new_group)=points_in_group(new_group)+1
            inside(point)=ins(p)
c            write(43,*)new_group,high_point,point,p
          end if
        end do
c     
        particle=hoc_particles(high_point)
c     
        if(particles_at_point(high_point)+highest_used_particle 
     >    .ge. particles_maxx) then
          memory_value=2
          print*,'exceeded particles in daughter ',highest_used_particle
          return
        end if
c     
        do while(particle .gt. 0)
          new_particle=new_particle+1
c     
          p_x=(pos_x(particle_pointer(particle))*grid_multiply-
     >      pos_point_x(point_tmp(1)))*d_point_inv
          p_y=(pos_y(particle_pointer(particle))*grid_multiply-
     >      pos_point_y(point_tmp(1)))*d_point_inv
          p_z=(pos_z(particle_pointer(particle))*grid_multiply-
     >      pos_point_z(point_tmp(1)))*d_point_inv
c     
          p=1+p_x+3*p_y+9*p_z
          n_p=point_tmp(p)
          next_particles(new_particle)=hoc_particles(n_p)
          hoc_particles(n_p)=new_particle
          particles_at_point(n_p)=particles_at_point(n_p)+1
          particle_pointer(new_particle)=particle_pointer(particle)
          particles_in_group(new_group)=
     >      particles_in_group(new_group)+1
c     
          particle=next_particles(particle)
        end do
c     
        highest_used_particle=new_particle
c
        high_point=next_high_points(high_point)
      end do
c     
      highest_used_point=point
c     
      call cpu_time(t_2)
c     
      delta=2**(level_max-level(new_group))
      high_point=hoc_high_points(high_group)
      do while(high_point .gt. 0)
c        print*,'neighbors 3 ',high_point
        call find_neighbors(high_point,adj,
     >    high_point_up_x,high_point_up_y,high_point_up_z,
     >    high_point_down_x,high_point_down_y,high_point_down_z,
     >    points_maxx)
c
        if(debug) call just_checking(adj,
     >    pos_point_x,pos_point_y,pos_point_z,
     >    hoc_daughter_points,next_daughter_points,delta,points_maxx)

c
        call look_where_erika(adj,high_point,
     >    high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
c
        call find_uppers_erika(adj,periodic,grid_multiply,
     >    high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    high_point,hoc_daughter_points,next_daughter_points,
     >    pos_point_x,pos_point_y,pos_point_z,delta,
     >    point_up_x,point_up_y,point_up_z,point_pointer,
     >    hoc_high,next_high_points,debug,high_points_maxx,points_maxx)
c
        call look_where_erika(adj,high_point,
     >    high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    decisions_xm,decisions_ym,decisions_zm,
     >    positions_xm,positions_ym,positions_zm)
c
        call find_downers_nina(adj,periodic,grid_multiply,
     >    high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    high_point,hoc_daughter_points,next_daughter_points,
     >    pos_point_x,pos_point_y,pos_point_z,delta,
     >    point_down_x,point_down_y,point_down_z,point_pointer,
     >    hoc_high,next_high_points,debug,high_points_maxx,points_maxx)
c     
        high_point=next_high_points(high_point)
      end do
c     
      call cpu_time(t_3)
c     
c      point=hoc_points(new_group)
c      do while(point .gt. 0)
c        points_in_group(new_group)=points_in_group(new_group)+1
c     
c        if(point_up_x(point) .gt. 0)
c     >    point_down_x(point_up_x(point))=point
c        if(point_up_y(point) .gt. 0)
c     >    point_down_y(point_up_y(point))=point
c        if(point_up_z(point) .gt. 0)
c     >    point_down_z(point_up_z(point))=point
c        point=next_points(point)
c      end do
c     
      call cpu_time(t_4)
c     
      total_points=0
      high_point=hoc_high_points(high_group)
      do while(high_point .gt. 0)
        point=hoc_daughter_points(high_point)
        do while(point .gt. 0)
          total_points=total_points+1
          if(point_pointer(point) .eq. -1) then
            point_pointer(point)=high_point
          else if(point_pointer(point) .eq. -3) then
            point_pointer(point)=point_up_x(high_point)
          else if(point_pointer(point) .eq. -7) then
            point_pointer(point)=point_up_y(high_point)
          else if(point_pointer(point) .eq. -9) then
            point_pointer(point)=point_up_x(point_up_y(high_point))
          else if(point_pointer(point) .eq. -19) then
            point_pointer(point)=point_up_z(high_point)
          else if(point_pointer(point) .eq. -21) then
            point_pointer(point)=point_up_x(point_up_z(high_point))
          else if(point_pointer(point) .eq. -25) then
            point_pointer(point)=point_up_y(point_up_z(high_point))
          else if(point_pointer(point) .eq. -27) then
            point_pointer(point)=
     >        point_up_y(point_up_x(point_up_z(high_point)))
          end if
          point=next_daughter_points(point)
        end do
        high_point=next_high_points(high_point)
      end do
c     
      if(mod(tweaks/2,2) .eq. 0 .and. debug) then
        point=hoc_points(new_group)
        do while(point .gt. 0)
          if(point_pointer(point) .gt. 0) then
            if(.not. inside(point_pointer(point))) then
              if(mother_group(new_group) .eq. 1 .and. 
     >          .not. periodic) then
                print*,'touching moat ',new_group,point,
     >            point_pointer(point)
              else
                print*,'buffer not working ',new_group,
     >            mother_group(new_group),level(new_group),
     >            point,point_pointer(point),
     >            pos_point_x(point),
     >            pos_point_y(point),
     >            pos_point_z(point)
              end if
            end if
          end if
          point=next_points(point)
        end do
      end if
c     
      call cpu_time(t_5)
      if(t_5-t_1 .gt. 0.1)
     >  write(39,39)new_group,total_points,t_5-t_1,t_2-t_1,t_3-t_2,
     >  t_4-t_3,t_5-t_4
 39   format(2i7,7(1pe12.4))
      return
      end
