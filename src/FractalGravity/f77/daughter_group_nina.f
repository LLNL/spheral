      subroutine daughter_group(new_group,hoc_points,
     >  next_points,hoc_high_points,next_high_points,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  hoc_particles,next_particles,pos_x,pos_y,
     >  pos_z,mother_group,high_group,inside,periodic,
     >  level,pos_point_x,pos_point_y,pos_point_z,
     >  point_pointer,level_max,it_is_high,
     >  particles_at_point,particles_in_group,
     >  points_in_group,particle_pointer,real_pointer,
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
      logical first_time,re_order,remove_duplicates
      logical decisions(27,27),it_is_high(-1:points_maxx)
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
      integer n_p,real_pointer(points_maxx)
      integer next_points(points_maxx),pos_point_x(points_maxx)
      integer pos_point_y(points_maxx),pos_point_z(points_maxx)
      integer hoc_particles(points_maxx),particles_at_point(points_maxx)
      integer point_pointer(points_maxx),particles_in_group(groups_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx),particle,new_particle
      integer next_particles(particles_maxx)
      integer daughter(high_points_max_group),n_high
      integer points_in_group(groups_maxx)
      integer grid_multiply,point_000,point_100,point_010,point_110
      integer point_001,point_101,point_011,point_111
      integer particle_pointer(particles_maxx)
      integer next_high_points(high_points_maxx)
      integer lowest_point_in_group(groups_maxx)
      integer highest_point_in_group(groups_maxx)
      integer memory_value
      integer adj(-1:1,-1:1,-1:1)
      integer positions(27,27),p,p_x,p_y,p_z,point_tmp(27)
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),d_point_inv
c     
      real t_1,t_2,t_3,t_4,t_5
c     
      data first_time/.true./
      data decisions/729*.false./
      data positions/729*-1000/
c     
      call cpu_time(t_1)
c     
      if(first_time) then
        first_time=.false.
        call make_decisions_erika(decisions,positions)
      end if
c     
c     write(59,*)'new group ',new_group
c     
      memory_value=0
      level(new_group)=level(mother_group(new_group))+1
      d_point=2**(level_max-level(new_group))
      d_point_inv=1.0/float(d_point)
      grid_multiply=grid_length*2**level_max
      delta=2**(level_max-level(new_group))
      new_particle=highest_used_particle
      point=highest_used_point
      hoc_points(new_group)=-1
      points_in_group(new_group)=0
      particles_in_group(new_group)=0
c     
      high_point=hoc_high_points(high_group)
      hoc_high=high_point
c     
      n_high=0
      do while(high_point .gt. 0)
c     
        n_high=n_high+1
c     
        if(point+27 .ge. points_maxx) then
          memory_value=1
          print*,'memory value= ',memory_value
          print*,'exceeded points in daughter ',point+27
          return
        end if
c     
        call find_neighbors(high_point,adj,it_is_high,
     >    point_up_x,point_up_y,point_up_z,
     >    point_down_x,point_down_y,point_down_z,
     >    points_maxx)
c     
        call find_go_ahead(go_ahead,ins,adj,decisions,positions,
     >    new_group)
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
            if(periodic) then
              pos_point_x(point)=mod(pos_point_x(point),grid_multiply)
              pos_point_y(point)=mod(pos_point_y(point),grid_multiply)
              pos_point_z(point)=mod(pos_point_z(point),grid_multiply)
            end if
            real_pointer(point)=p
            hoc_particles(point)=-1
            particles_at_point(point)=0
            point_pointer(point)=-p
            if(p .eq. 1) then
              daughter(n_high)=point
            end if
            point_tmp(p)=point
            points_in_group(new_group)=points_in_group(new_group)+1
            inside(point)=ins(p)
          end if
        end do
c     
        particle=hoc_particles(high_point)
c     
        if(particles_at_point(high_point)+highest_used_particle 
     >    .ge. particles_maxx) then
          memory_value=2
          print*,'exceeded particles in daughter ',
     >      highest_used_particle
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
        high_point=next_high_points(high_point)
c     
      end do
c     
      highest_used_point=point
c     
      call cpu_time(t_2)
c     
      re_order=.true.
      remove_duplicates=.true.
c     
      call my_list_sort_3(hoc_points(new_group),next_points,
     >  pos_point_x,pos_point_y,pos_point_z,point_up_x,point_pointer,
     >  delta,periodic,grid_multiply,remove_duplicates,re_order,
     >  points_maxx)
c     
      re_order=.false.
      remove_duplicates=.false.
c     
      call my_list_sort_3(hoc_points(new_group),next_points,
     >  pos_point_z,pos_point_x,pos_point_y,point_up_z,point_pointer,
     >  delta,periodic,grid_multiply,remove_duplicates,re_order,
     >  points_maxx)
c     
      re_order=.false.
      remove_duplicates=.false.
c     
      call my_list_sort_3(hoc_points(new_group),next_points,
     >  pos_point_y,pos_point_z,pos_point_x,point_up_y,point_pointer,
     >  delta,periodic,grid_multiply,remove_duplicates,re_order,
     >  points_maxx)
c     
      call cpu_time(t_3)
c     
      point=hoc_points(new_group)
      do while(point .gt. 0)
        point_down_x(point)=-1
        point_down_y(point)=-1
        point_down_z(point)=-1
        point=next_points(point)
      end do
c     
      point=hoc_points(new_group)
      do while(point .gt. 0)
c     
        if(point_up_x(point) .gt. 0)
     >    point_down_x(point_up_x(point))=point
        if(point_up_y(point) .gt. 0)
     >    point_down_y(point_up_y(point))=point
        if(point_up_z(point) .gt. 0)
     >    point_down_z(point_up_z(point))=point
        point=next_points(point)
      end do
c     
      call test_up_down(hoc_points(new_group),new_group,
     >  level(new_group),next_points,
     >  point_pointer,point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  pos_point_x,pos_point_y,pos_point_z,points_maxx)
c     
      call cpu_time(t_4)
c     
      high_point=hoc_high_points(high_group)
      n_high=0
      do while(high_point .gt. 0)
        n_high=n_high+1
        point_000=daughter(n_high)
c     
        point_100=point_up_x(point_up_x(point_000))
        point_010=point_up_y(point_up_y(point_000))
        point_110=point_up_y(point_up_y(point_100))
        point_001=point_up_z(point_up_z(point_000))
        point_101=point_up_x(point_up_x(point_001))
        point_011=point_up_y(point_up_y(point_001))
        point_111=point_up_y(point_up_y(point_101))
c     
        point_pointer(point_000)=high_point
        point_pointer(point_100)=point_up_x(high_point)
        point_pointer(point_010)=point_up_y(high_point)
        point_pointer(point_110)=point_up_y(point_pointer(point_100))
        point_pointer(point_001)=point_up_z(high_point)
        point_pointer(point_101)=point_up_x(point_pointer(point_001))
        point_pointer(point_011)=point_up_y(point_pointer(point_001))
        point_pointer(point_111)=point_up_y(point_pointer(point_101))
c     
        high_point=next_high_points(high_point)
      end do
c     
      call cpu_time(t_5)
      if(t_5-t_1 .gt. 0.1)
     >  write(39,39)new_group,points_in_group(new_group),
     >  t_5-t_1,t_2-t_1,t_3-t_2,t_4-t_3,t_5-t_4
 39   format(2i9,7(1pe12.4))
      return
      end
