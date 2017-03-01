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
c
      logical do_it(27)
      logical periodic,inside(points_maxx),debug,what_1,what_2,what_3
      logical go_ahead,empty(-1:1,-1:1,-1:1),test_it,first_time
      logical decisions_x(27,27),decisions_y(27,27),decisions_z(27,27)
      logical decisions(27,27)
c     
      include 'maxx.inc'
c     
      integer mother_group(groups_maxx),new_group,pos
      integer d_point,level_max,level(groups_maxx),grid_length
      integer highest_used_particle,point,highest_used_point
      integer hoc_points(groups_maxx),high_point,delta
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
      integer position(19)
      integer n,i,j
      integer adj(-1:1,-1:1,-1:1),total_points,list_high
      integer positions(27,27)
      integer positions_x(27,27),positions_y(27,27),positions_z(27,27)
c
      integer high_x(maxx),high_y(maxx),high_z(maxx)
      integer pos_high_x(maxx),pos_high_y(maxx),pos_high_z(maxx)
c
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),d_point_inv
c     
      real t_1,t_2,t_3,t_4,t_5
c     
      test_it(i,j)=i .lt. 0 .or. i .gt. j
c     
      data first_time/.true./
      data delta_x/2,2,2,2,0,1,0,1,2,2,0,1,0,1,2,2,0,1,2/
      data delta_y/0,1,0,1,2,2,2,2,2,2,0,0,1,1,0,1,2,2,2/
      data delta_z/0,0,1,1,0,0,1,1,0,1,2,2,2,2,2,2,2,2,2/
      data position/3,6,12,15,7,8,16,17,9,18,19,20,22,23,21,24,25,26,27/
      data decisions_x/729*.false./
      data decisions_y/729*.false./
      data decisions_z/729*.false./
      data positions_x/729* -1000/
      data positions_y/729* -1000/
      data positions_z/729* -1000/
      data do_it/
     >  .true.,.true.,.false.,.true.,.true.,.false.,
     >  .false.,.false.,.false.,.true.,.true.,.false.,
     >  .true.,.true.,.false.,.false.,.false.,.false.,
     >  .false.,.false.,.false.,.false.,.false.,.false.,
     >  .false.,.false.,.false./
c     
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
      if(first_time) then
        first_time=.false.
        call make_decisions_0(decisions,positions)
        call make_decisions(decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
      end if
c     
      list_high=0
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
        list_high=list_high+1
c     
        call find_empties(high_point,empty,
     >    high_point_up_x,high_point_up_y,high_point_up_z,
     >    high_point_down_x,high_point_down_y,high_point_down_z,
     >    points_maxx)
c     
        call find_neighbors(high_point,adj,
     >    high_point_up_x,high_point_up_y,high_point_up_z,
     >    high_point_down_x,high_point_down_y,high_point_down_z,
     >    points_maxx)
c
        write(42,*)'here ',new_group,high_point
        call look_where(adj,list_high,do_it,
     >    high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    decisions_x,decisions_y,decisions_z,
     >    positions_x,positions_y,positions_z)
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
          go_ahead=do_it(position(n))
c          if(n .eq. 1) then
c            pos=3
c            go_ahead=empty(1,0,0)
c          else if(n .eq. 2) then
c            pos=6
c            go_ahead=empty(1,0,0)
c          else if(n .eq. 3) then
c            pos=12
c            go_ahead=empty(1,0,0)
c          else if(n .eq. 4) then
c            pos=15
c            go_ahead=empty(1,0,0)
c          else if(n .eq. 5) then
c            pos=7
c            go_ahead=empty(0,1,0) .and. empty(-1,1,0)
c          else if(n .eq. 6) then
c            pos=8
c            go_ahead=empty(0,1,0)
c          else if(n .eq. 7) then
c            pos=16
c            go_ahead=empty(0,1,0) .and. empty(-1,1,0)
c          else if(n .eq. 8) then
c            pos=17
c            go_ahead=empty(0,1,0)
c          else if(n .eq. 9) then
c            pos=9
c            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
c     >        empty(1,1,0)
c          else if(n .eq. 10) then
c            pos=18
c            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
c     >        empty(1,1,0)
c          else if(n .eq. 11) then
c            pos=19
c            go_ahead=empty(0,0,1) .and. empty(-1,0,1) .and.
c     >        empty(0,-1,1) .and. empty(-1,-1,1)
c          else if(n .eq. 12) then
c            pos=20
c            go_ahead=empty(0,0,1) .and. empty(0,-1,1)
c          else if(n .eq. 13) then
c            pos=22
c            go_ahead=empty(0,0,1) .and. empty(-1,0,1)
c          else if(n .eq. 14) then
c            pos=23
c            go_ahead=empty(0,0,1)
c          else if(n .eq. 15) then
c            pos=21
c            go_ahead=empty(1,0,0) .and. 
c     >        empty(0,0,1) .and. empty(1,0,1) .and.
c     >        empty(0,-1,1) .and. empty(1,-1,1)
c          else if(n .eq. 16) then
c            pos=24
c            go_ahead=empty(1,0,0) .and. 
c     >        empty(0,0,1) .and. empty(1,0,1)
c          else if(n .eq. 17) then
c            pos=25
c            go_ahead=empty(0,1,0) .and. empty(-1,1,0) .and.
c     >        empty(0,0,1) .and. empty(0,1,1) .and.
c     >        empty(-1,0,1) .and. empty(-1,1,1)
c          else if(n .eq. 18) then
c            pos=26
c            go_ahead=empty(0,1,0) .and.
c     >        empty(0,0,1) .and. empty(0,1,1)
c          else if(n .eq. 19) then
c            pos=27
c            go_ahead=empty(1,0,0) .and. empty(0,1,0) .and.
c     >        empty(1,1,0) .and.
c     >        empty(0,0,1) .and. empty(1,0,1) .and.
c     >        empty(0,1,1) .and. empty(1,1,1)
c          end if
c     
          if(go_ahead) then
            point=point+1
            point_pointer(point)=-position(n)
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
      delta=2**(level_max-level(new_group))
      list_high=0.0
      high_point=hoc_high_points(high_group)
      do while(high_point .gt. 0)
        list_high=list_high+1
        call find_uppers_hope(list_high,high_x,high_y,high_z,
     >    pos_high_x,pos_high_y,pos_high_z,
     >    high_point,hoc_daughter_points,next_daughter_points,
     >    pos_point_x,pos_point_y,pos_point_z,delta,
     >    point_up_x,point_up_y,point_up_z,point_pointer,points_maxx)
c
        high_point=next_high_points(high_point)
      end do
c     
      call cpu_time(t_3)
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
