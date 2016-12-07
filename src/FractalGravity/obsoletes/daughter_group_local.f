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
      integer groups_maxx,points_maxx,points_maxx_group,high_points_maxx
      integer particles_maxx,particles_real,tweaks
      logical periodic,inside(points_maxx),debug,what_1,what_2,what_3
      logical go_ahead,empty(-1:1,-1:1,-1:1),test_it
      logical check_remove
c     
      include 'maxx.inc'
c     
      integer mother_group(groups_maxx),new_group,pos
      integer d_point,level_max,level(groups_maxx),grid_length
      integer highest_used_particle,point,highest_used_point
      integer hoc_points(groups_maxx),high_point
      integer hoc_high_points(groups_maxx),high_group
      integer point_0,point_down_x(points_maxx)
      integer point_down_y(points_maxx),point_down_z(points_maxx)
      integer n_p,n_p_x,n_p_y,n_p_z
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
      integer delta_x(19),delta_y(19),delta_z(19),memory_value
      integer n,i,j,hoc_possibles(maxx),next_possibles(maxx)
      integer hoc_backwards,next_backwards(maxx)
      integer adj(-1:1,-1:1,-1:1)
c     
      integer total_points
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),d_point_inv
c     
      real t_1,t_2,t_3,t_4,t_5,t_6,t_7
c     
      test_it(i,j)=i .lt. 0 .or. i .gt. j
c     
      data delta_x/2,2,2,2,0,1,0,1,2,2,0,1,0,1,2,2,0,1,2/
      data delta_y/0,1,0,1,2,2,2,2,2,2,0,0,1,1,0,1,2,2,2/
      data delta_z/0,0,1,1,0,0,1,1,0,1,2,2,2,2,2,2,2,2,2/
c     
c***  call FIND_EMPTIES
c***  find empty neighbor virtual high cells
c***  call FIND_UPPERS
c***  find point_up_x etc

c     
      call cpu_time(t_1)
c     
c     print*,' I am here in daughter group'
c     
      memory_value=0
      go_ahead=.false.
      level(new_group)=level(mother_group(new_group))+1
      d_point=2**(level_max-level(new_group))
      d_point_inv=1.0/float(d_point)
      grid_multiply=grid_length*2**level_max
      new_particle=highest_used_particle
      point=highest_used_point
      hoc_points(new_group)=-1
      pos=-1000
      points_in_group(new_group)=0
      particles_in_group(new_group)=0
c     
      call find_possibles(high_group,new_group,
     >  hoc_high_points,next_high_points,
     >  pos_point_x,grid_length,level_max,level,hoc_possibles,
     >  next_possibles,groups_maxx,points_maxx)
c     
      high_point=hoc_high_points(high_group)
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
        hoc_daughter_points(high_point)=-1
        point_0=point
c     
c     call find_real_empties(high_group,high_point,
c     >    hoc_high_points,next_high_points,periodic,
c     >    grid_multiply,d_point,pos_point_x,pos_point_y,pos_point_z,
c     >    empty,groups_maxx,points_maxx)
c     
        call find_empties(high_point,empty,
     >    high_point_up_x,high_point_up_y,high_point_up_z,
     >    high_point_down_x,high_point_down_y,high_point_down_z)
c     
c***  generate the points at positions 1,2,4,5,10,11,13,14
        do n_p_z=0,1
          do n_p_y=0,1
            do n_p_x=0,1
              point=point+1
              next_points(point)=hoc_points(new_group)
              hoc_points(new_group)=point
              pos_point_x(point)=pos_point_x(high_point)+
     >          n_p_x*d_point
              pos_point_y(point)=pos_point_y(high_point)+
     >          n_p_y*d_point
              pos_point_z(point)=pos_point_z(high_point)+
     >          n_p_z*d_point
              hoc_particles(point)=-1
              particles_at_point(point)=0
              next_daughter_points(point)=
     >          hoc_daughter_points(high_point)
              hoc_daughter_points(high_point)=point
            end do
          end do
        end do
c     
        point_pointer(point_0+1)=-1
        point_pointer(point_0+2)=-2
        point_pointer(point_0+3)=-4
        point_pointer(point_0+4)=-5
c     
        point_pointer(point_0+5)=-10
        point_pointer(point_0+6)=-11
        point_pointer(point_0+7)=-13
        point_pointer(point_0+8)=-14
c     
        inside(point_0+1)=.not.(
     >    empty(-1,0,0) .or. empty(0,-1,0) .or. empty(-1,-1,0) .or.
     >    empty(-1,0,-1) .or. empty(0,-1,-1) .or. empty(-1,-1,-1)
     >    .or. empty(0,0,-1))
        inside(point_0+2)=.not.(
     >    empty(0,0,-1) .or. empty(0,-1,0) .or. empty(0,-1,-1))
        inside(point_0+3)=.not.(
     >    empty(0,0,-1) .or. empty(-1,0,0) .or. empty(-1,0,-1))
        inside(point_0+4)=.not. empty(0,0,-1)
        inside(point_0+5)=.not.(
     >    empty(-1,0,0) .or. empty(0,-1,0) .or. empty(-1,-1,0))
        inside(point_0+6)=.not. empty(0,-1,0)
        inside(point_0+7)=.not. empty(-1,0,0)
        inside(point_0+8)=.true.
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
          n_p_x=(pos_x(particle_pointer(particle))*grid_multiply-
     >      pos_point_x(point_0+1))*d_point_inv
          n_p_y=(pos_y(particle_pointer(particle))*grid_multiply-
     >      pos_point_y(point_0+1))*d_point_inv
          n_p_z=(pos_z(particle_pointer(particle))*grid_multiply-
     >      pos_point_z(point_0+1))*d_point_inv
c     
          n_p=1+n_p_x+2*n_p_y+4*n_p_z+point_0
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
        do n=1,19
          if(n .eq. 1) then
            pos=3
            go_ahead=empty(1,0,0)
          else if(n .eq. 2) then
            pos=6
            go_ahead=empty(1,0,0)
          else if(n .eq. 3) then
            pos=12
            go_ahead=empty(1,0,0)
          else if(n .eq. 4) then
            pos=15
            go_ahead=empty(1,0,0)
          else if(n .eq. 5) then
            pos=7
            go_ahead=empty(0,1,0) .and. empty(-1,1,0)
          else if(n .eq. 6) then
            pos=8
            go_ahead=empty(0,1,0)
          else if(n .eq. 7) then
            pos=16
            go_ahead=empty(0,1,0) .and. empty(-1,1,0)
          else if(n .eq. 8) then
            pos=17
            go_ahead=empty(0,1,0)
          else if(n .eq. 9) then
            pos=9
            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
     >        empty(1,1,0)
          else if(n .eq. 10) then
            pos=18
            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
     >        empty(1,1,0)
          else if(n .eq. 11) then
            pos=19
            go_ahead=empty(0,0,1) .and. empty(-1,0,1) .and.
     >        empty(0,-1,1) .and. empty(-1,-1,1)
          else if(n .eq. 12) then
            pos=20
            go_ahead=empty(0,0,1) .and. empty(0,-1,1)
          else if(n .eq. 13) then
            pos=22
            go_ahead=empty(0,0,1) .and. empty(-1,0,1)
          else if(n .eq. 14) then
            pos=23
            go_ahead=empty(0,0,1)
          else if(n .eq. 15) then
            pos=21
            go_ahead=empty(1,0,0) .and. 
     >        empty(0,0,1) .and. empty(1,0,1) .and.
     >        empty(0,-1,1) .and. empty(1,-1,1)
          else if(n .eq. 16) then
            pos=24
            go_ahead=empty(1,0,0) .and. 
     >        empty(0,0,1) .and. empty(1,0,1)
          else if(n .eq. 17) then
            pos=25
            go_ahead=empty(0,1,0) .and. empty(-1,1,0) .and.
     >        empty(0,0,1) .and. empty(0,1,1) .and.
     >        empty(-1,0,1) .and. empty(-1,1,1)
          else if(n .eq. 18) then
            pos=26
            go_ahead=empty(0,1,0) .and.
     >        empty(0,0,1) .and. empty(0,1,1)
          else if(n .eq. 19) then
            pos=27
            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
     >        empty(1,1,0) .and.
     >        empty(0,0,1) .and. empty(1,0,1) .and.
     >        empty(0,1,1) .and. empty(1,1,1)
          end if
c     
          if(go_ahead) then
            point=point+1
            point_pointer(point)=-pos
            next_points(point)=hoc_points(new_group)
            hoc_points(new_group)=point
            hoc_particles(point)=-1
            particles_at_point(point)=-1
            inside(point)=.false.
            next_daughter_points(point)=
     >        hoc_daughter_points(high_point)
            hoc_daughter_points(high_point)=point
            pos_point_x(point)=pos_point_x(high_point)+
     >        delta_x(n)*d_point
            pos_point_y(point)=pos_point_y(high_point)+
     >        delta_y(n)*d_point
            pos_point_z(point)=pos_point_z(high_point)+
     >        delta_z(n)*d_point
            if(periodic .and. debug) then
              what_1=test_it(pos_point_x(point),grid_multiply)
              what_2=test_it(pos_point_y(point),grid_multiply)
              what_3=test_it(pos_point_z(point),grid_multiply)
              if(what_1 .or. what_2 .or. what_3) print*,'outsider'
            end if
          end if
        end do
        high_point=next_high_points(high_point)
      end do
c     
      highest_used_point=point
c     
      call cpu_time(t_2)
c     
      check_remove=next_high_points(hoc_high_points(high_group)) .gt. 0
c     
c      check_remove=.false.
      if(check_remove) then
c     
        call backwards(new_group,hoc_points,next_points,
     >    hoc_backwards,next_backwards)
        call cpu_time(t_3)
c     if(debug) then
c     point=hoc_points(new_group)
c     do while(point .gt. 0)
c     write(81,*)'first  ',new_group,point,
c     >        pos_point_x(point),pos_point_y(point),pos_point_z(point)
c     point=next_points(point)
c     end do
c     end if
c     
        high_point=hoc_high_points(high_group)
        do while(high_point .gt. 0)
          call remove_duplicates(high_group,high_point,
     >      hoc_high_points,next_high_points,hoc_points,next_points,
     >      hoc_daughter_points,next_daughter_points,new_group,
     >      periodic,grid_multiply,d_point,hoc_backwards,next_backwards,
     >      pos_point_x,pos_point_y,pos_point_z,point_pointer,
     >      hoc_possibles,next_possibles,
     >      inside,groups_maxx,points_maxx)
          
          high_point=next_high_points(high_point)
        end do
        
        call cpu_time(t_4)
      else
        t_3=t_2
        t_4=t_3
      end if
c     
c     if(debug) then
c     point=hoc_points(new_group)
c     do while(point .gt. 0)
c     write(81,*)'second ',new_group,point,
c     >      pos_point_x(point),pos_point_y(point),pos_point_z(point)
c     point=next_points(point)
c     end do
c     end if
c     
      high_point=hoc_high_points(high_group)
      do while(high_point .gt. 0)
c****
c        call find_neighbors(high_point,adj,
c     >    high_point_up_x,high_point_up_y,high_point_up_z,
c     >    high_point_down_x,high_point_down_y,high_point_down_z,
c     >    points_maxx)
c***
c        call find_uppers(high_point,hoc_daughter_points,
c     >    high_point_up_x,high_point_up_y,high_point_up_z,
c     >    next_daughter_points,point_pointer,adj,points_maxx)
c***
        call find_uppers(high_group,high_point,hoc_high_points,
     >    next_high_points,hoc_daughter_points,next_daughter_points,
     >    periodic,grid_multiply,d_point,
     >    point_up_x,point_up_y,point_up_z,
     >    pos_point_x,pos_point_y,pos_point_z,point_pointer,
     >    groups_maxx,points_maxx)
c     
        high_point=next_high_points(high_point)
      end do
c     
      call cpu_time(t_5)
c     
      point=hoc_points(new_group)
      do while(point .gt. 0)
        points_in_group(new_group)=points_in_group(new_group)+1
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
      call cpu_time(t_6)
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
      call cpu_time(t_7)
      if(t_7-t_1 .gt. 1.0)
     >  write(39,39)new_group,total_points,t_7-t_1,t_2-t_1,t_3-t_2,
     >  t_4-t_3,t_5-t_4,t_6-t_5,t_7-t_6
 39   format(2i7,7(1pe12.4))
      return
      end
