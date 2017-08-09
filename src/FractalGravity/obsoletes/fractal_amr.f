      subroutine fractal_amr(number_particles,grid_length,periodic,
     >     minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >     potential,force_x,force_y,force_z,
     >     moat_0,random_offset,debug,tweaks,memory_value)
c
c*** INPUT IS UNCHANGED EXCEPT PERIODIC DATA IS WRAPPED AROUND
c*** POSITIONS MAY BE OFFSET BUT THEY ARE PUT BACK BEFORE RETURNING
      implicit none
c
      integer groups_max,points_max,particles_max,number_particles
      integer points_max_group,high_points_max,iso_maxx,maxx
c
      include 'maxx.inc'
      include 'maxes.inc'
c
      logical inside(points_max),periodic,work_to_be_done,debug
      logical done_group(groups_max),chain_finished(groups_max)
      logical chain_started(groups_max),ok_start_chain(groups_max)
      real pos_x(number_particles),pos_y(number_particles)
      real pos_z(number_particles),particle_mass(number_particles)
      real potential_point(points_max)
      real potential(number_particles),density(points_max)
      real force_point_x(points_max),force_point_y(points_max)
      real force_point_z(points_max)
      real force_x(number_particles),force_y(number_particles)
      real force_z(number_particles)
c
      integer random_offset,number_chains,memory_value,tweaks
      integer grid_length,level_max,group,group_loop
      integer next_points(points_max),hoc_points(groups_max)
      integer hoc_particles(points_max),next_particles(particles_max)
      integer pos_point_x(points_max),pos_point_y(points_max)
      integer pos_point_z(points_max)
      integer point_up_x(points_max),point_up_y(points_max)
      integer point_up_z(points_max)
      integer point_down_x(points_max),point_down_y(points_max)
      integer point_down_z(points_max),point_pointer(points_max)
      integer high_point_up_x(high_points_max)
      integer high_point_up_y(high_points_max)
      integer high_point_up_z(high_points_max)
      integer high_point_down_x(high_points_max)
      integer high_point_down_y(high_points_max)
      integer high_point_down_z(high_points_max)
      integer particles_at_point(points_max)
      integer mother_of_high_group(groups_max)
      integer highest_used_point,highest_used_high_point
      integer highest_used_particle,highest_used_group,new_group
      integer particle_pointer(particles_max)
      integer minimum_number,seq_to_list(points_max_group),high_pairs
      integer list_high_1(points_max_group)
      integer list_high_2(points_max_group)
      integer group_tmp(points_max_group)
      integer next_groups(groups_max),hoc_groups
      integer level(groups_max),level_high(groups_max)
      integer number_high_groups(groups_max)
      integer highest_used_high_group,number_high_points
      integer mother_group_current,hoc_high_groups(groups_max)
      integer next_high_groups(groups_max),mother_group(groups_max)
      integer hoc_high_points(groups_max)
      integer next_high_points(high_points_max)
      integer hoc_groups_chain(groups_max),next_groups_chain(groups_max)
      integer highest_level_group(particles_max),grid_multiply
      integer groups_maxx,points_maxx,particles_maxx,particles_real
      integer high_points_maxx,points_maxx_group,moat,moat_0,n
      integer high_group,chain
c
c***  groups_max(I) maximum number of groups
c***  points_max(I) maximum number of points
c***  high_points_max(I) maximum number of high points
c***     currently points_max=high_points_max
c***  points_max_group(I) maximum number of points in a single group
c***  particles_max(I) maximum number of particles
c
c***  random_offset(I) random number seed
c***                   if .ne. 0 offset particles to avoid interference with lattice
c***                   if = 0 do nothing 
c***  inside(point;L) is this an inside point?
c***  periodic(L)  periodic boundary condition or isolate boundary conditions
c***  work_to_be_done (L) are there still groups possibly not done yet
c***  pos_x(particle;R) x position of particles 0< pos_x < 1
c***  pos_y(particle;R) y position of particles 0< pos_y < 1
c***  potential_point(point;R) potential at a point
c***  potential(particle;R) potential at a particle
c***  density(point;R) density at a point
c***  force_point_x(point;R) x-force at a point
c***  force_point_y(point;R) y-force at a point
c***  force_x(particle;R) x-force at a particle
c***  force_y(particle;R) y-force at a particle
c***  force_z(particle;R) z-force at a particle
c***  grid_multiply(I) effective number of points across
c***  grid_length(I) number of grid cells across in each dimension
c***  level_max(I) maximum level of a group
c***  group(I) current group I am working on
c***  group_loop(I) loop variable for groups
c***  next_points(point;I) next point in the point link list
c***  hoc_points(group;I) hoc of points in group
c***  hoc_particles(point;I) hoc of particles at point
c***  next_particles(particle;I) next particle in particle link list
c***  pos_point_x(point;I) x-position of point 
c***             0 <= pos_point_x < grid_length * 2**level_max
c***  pos_point_y(point;I) y-position of point 
c***             0 <= pos_point_y < grid_length * 2**level_max
c***  pos_point_z(point;I) z-position of point 
c***             0 <= pos_point_z < grid_length * 2**level_max
c***  point_up_x(point;I) neighbor point up in x (if < 0 no neigh point)
c***  point_up_y(point;I) neighbor point up in y (if < 0 no neigh point)
c***  point_up_z(point;I) neighbor point up in z (if < 0 no neigh point)
c***  point_down_x(point;I) neighbor point down in x (if < 0 no neigh point)
c***  point_down_y(point;I) neighbor point down in y (if < 0 no neigh point)
c***  point_down_z(point;I) neighbor point down in z (if < 0 no neigh point)
c***  high_point_up_x(high_point;I) 
c***             neighbor high point up in x (if < 0 no high neigh point)
c***  high_point_up_y(high_point;I) 
c***             neighbor high point up in y (if < 0 no high neigh point)
c***  high_point_up_z(high_point;I) 
c***             neighbor high point up in z (if < 0 no high neigh point)
c***  high_point_down_x(high_point;I) 
c***             neighbor high point down in x (if < 0 no high neigh point)
c***  high_point_down_y(high_point;I) 
c***             neighbor high point down in y (if < 0 no high neigh point)
c***  high_point_down_z(high_point;I) 
c***             neighbor high point down in z (if < 0 no high neigh point)
c***  particles_a_point(point;I) number of particles at point
c***  mother_of_high_group(high_group;I) mother group of high group
c***  highest_used_point(I) highest_used_point+1 is first new point available
c***  highest_used_high_point(I) 
c***             highest_used_high_point+1 is first new point available
c***  highest_used_particle(I)  
c***             highest_used_particle+1 first new virtual particle available
c***  highest_used_group(I)  
c***             highest_used_group+1 first new group available
c***  new group(I)   new daughter group just generated
c***  particle_pointer(particle;I) points to original particle in group 1
c***  minimum_number(I) minimum number of particles to make a high point 
c***  seq_to_list(high_point,I) pointer from seq list 
c***                 of high points to link list of points
c***  high_pairs)(I) number of pairs of high points output by find_high_pairs
c***  list_high_1(high_point,I) 1st list of high points in high pairs
c***  list_high_2(high_point,I) 2nd list of high points in high pairs
c***  group_tmp(point;I) temp group number output by equivalence_class
c***  next_groups(group;I) next group in group link list
c***  hoc_groups(I)       hoc of group link list
c***  level(group;I)  level of group
c***  level_high(high_group;I)  level of high group
c***  number_high_groups(I)  total number of high groups
c***  highest_used_high_group(I)  guess!!
c***  number_high_points(I) total number of high points
c***  mother_group_current   guess!!
c***  hoc_high_groups(group;I) hoc of link list of high groups in group
c***  next_high_groups(high_group;I) next high group in link list
c***  mother_group(high_group;I) mother group of high group
c***  hoc_high_points(group;I) hoc of link list of high points in group
c***  next_high_points(point;I) next high point in link list 
c***                  of high points in group
c***  highest_level_group(particle,I) highest level group the 
c***                     particle belongs to
c***  groups_maxx(I)    =groups_max   used for passing arrays
c***  points_maxx(I)    =points_max   used for passing arrays
c***  points_maxx_group(I)    =points_max_group   used for passing arrays
c***  particles_maxx(I) =particles_max   used for passing arrays
c***  particles_real(I) =number_particles   used for passing arrays
c***  moat_0(I) input buffer in isolated calculation
c***  moat(I)=min(moat_0,1)  buffer in isolated calculation
c***  n(I)        dummy variable
c***  memory_value(I)  reurned by daughter_group
c***                  if =0 no memory problem
c***                  if =1 need more total points
c***                  if =2 need more particles
c***                  if =3 need more points per group
c***                  if =4 need more groups
c***                  anything else die
c
c
c
c
c*** Mr. Garibaldi to hapless Engineer: 
c*** It is a prototype, it does not have to be perfect,
c*** It just has to work.
c
c*** Sheridan to Mr. Garibaldi:
c*** I am not interested in your problems, 
c*** I am only interested in your solutions.
c
c*** Its fundamental design flaws are obscured
c*** by its superficial design flaws
c
c*** This is a prototype of fractal amr, it is not supposed to be pretty
c*** or fast, it is just supposed to work.
c*** This is a single node version
c
c*** It is simple to transform into a multi-node code
c*** run TREESTART first
c*** Build the tree with calls to FIND_HIGH_PAIRS, EQUIVALENCE_CLASS,
c*** HIGH_GROUPS and DAUGHTER_GROUP. Finish building the tree with calls to
c*** CLEAN_GROUPS_LIST and FIND_HIGHEST_LEVEL_GROUP
c*** While the tree is being built, call ASSIGN_DENSITY for a single group
c*** then call PERIODIC_SOLVER or ISOLATED_SOLVER.
c*** When this is all done the real work starts.
c*** call ASSIGN_DENSITY for all the groups,except #1, 
c*** these calls are independent so they can be farmed out to separate nodes
c*** For all groups at level 1 we independently call 
c*** POTENTIAL_START followed by POISSON_SOLVER and FORCE_AT_POINT
C*** then for level 2 and then for level 3 until we reach the maximum level.
c*** All these calls at a certain level are independent so they can be 
c*** farmed out to different nodes. The only restriction is the levels 
c*** must be done in ascending order. A more clever way is to say that 
c*** any group can start as soon as its mother group is finished.
c*** Once all groups are finished we can call
c*** FORCE_AT_PARTICLE 
c
c***  Input:
c***  number_particles (I)(Got that?)
c***  grid_length      (I)(number of cells in each dimension)
c***  periodic         (L) T/F is periodic/isolated boundary conditions
c***  minimum_number   (I) minimum number of particles in a cell so that 
c***                       it is a high density cell
c***  level_max        (I) maximum level of groups
c***  pos_x            (R) x-positions for particles [0,1]
c***  pos_y            (R) y-positions for particles [0,1]
c***  pos_z            (R) z-positions for particles [0,1]
c***                       if periodic, particles will be wrapped
c***                       if isolated, outside particles will be ignored
c***  output:
c***  potential        (R) potential at particle
c***  force_x          (R) x-force at particles
c***  force_y          (R) y-force at particles
c***  force_z          (R) z-force at particles
c***  moat_0           (I) size of moat for isolated boundary conditions
c***                       moat=min(moat_0,1). Not used in periodic calc
c***  random_offset    (I) if =0 no particle offset
c***                       if != seed for random offset
c***  debug            (L) T:debug run, F:no debug
c***  
c
c**** memory requirements (worst case scenario)
c***  total_particles=number_particles*(level_max+1)
c***  total points=grid_length**3+
c***  27/8*number_particles*
c***  (8/minimum_number-(8/minimum_number)**(level_max+1))/(1-8/minimum_number)
c***  if(minimum_number=8)
c***  total points=grid_length**3+27/8*number_particles*level_max
c
c**** Memory_worst=
c***          (L) 1 x total points
c***          (I) 22 x total points + 4 x points in biggest group
c***          (R) 5 x total points
c***          (I) 3 x total particles 
c***          (R) 8 x number_particles
c***          (I) 15 x number of groups
c***          (I,R,L) < 1000 eveything else
c
      call timing(-1,30)
c
      memory_value=0
      grid_multiply=grid_length*2**level_max
      groups_maxx=groups_max
      points_maxx=points_max
      high_points_maxx=high_points_max
      points_maxx_group=points_max_group
      particles_maxx=particles_max
      particles_real=number_particles
c

      moat=min(moat_0,1)
c
      if(periodic) then
         if((grid_length+1)**3 .gt. 
     >        min(points_max,points_max_group)) then
            memory_value=1
            print*,'memory value= ',memory_value
            return
         end if
      else
         if((2*grid_length+2)**3 .gt. 
     >        iso_maxx) then
            memory_value=1
            print*,'memory value= ',memory_value
            return
         end if
      end if
c
      if(number_particles .gt. particles_max) then
         memory_value=2
         print*,'memory value= ',memory_value
         return
      end if
c
      if(debug) write(31,*)'maxxes= ',
     >     particles_maxx,points_maxx,groups_maxx
c
      if(debug)write(31,*)'offset it'
      call timing(-1,1)
      call offset(1,random_offset,number_particles,
     >     grid_multiply,pos_x,pos_y,pos_z,tweaks)
      call timing(1,1)
c
      if(debug)write(31,*)'be negative'
      call timing(-1,2)
      call negatives_groups(groups_maxx,hoc_points,mother_of_high_group,
     >     number_high_groups,next_groups,level,level_high,
     >     hoc_high_groups,next_high_groups,mother_group,tweaks)
c
      call negatives_points(points_maxx,high_points_maxx,
     >     points_maxx_group,next_points,hoc_particles,pos_point_x,
     >     pos_point_y,pos_point_z,point_up_x,point_down_x,point_up_y,
     >     point_down_y,point_up_z,point_down_z,inside,
     >     particles_at_point,
     >     potential_point,force_point_x,force_point_y,force_point_z,
     >     high_point_up_x,high_point_down_x,
     >     high_point_up_y,high_point_down_y,
     >     high_point_up_z,high_point_down_z,seq_to_list,
     >     list_high_1,list_high_2,group_tmp,next_high_points,tweaks)
c
      call negatives_particles(particles_maxx,particles_real,
     >     next_particles,
     >     particle_pointer,highest_level_group,
     >     potential,force_x,force_y,force_z,tweaks)
      call timing(1,2)
c
      group=1
      next_groups(group)=-1
      hoc_groups=group
      level(group)=0
      number_high_groups(group)=-1
      highest_used_group=1
      highest_used_high_group=0
      highest_used_particle=0
      highest_used_point=0
      highest_used_high_point=0
c
      if(debug)print*,'treestart'
c
      call timing(-1,3)
      call tree_start(grid_length,level_max,group,groups_maxx,
     >     points_maxx,particles_maxx,periodic,hoc_points,
     >     next_points,highest_used_point,highest_used_particle,
     >     point_up_x,point_up_y,point_up_z,point_pointer,
     >     hoc_particles,next_particles,particle_pointer,
     >     point_down_x,point_down_y,point_down_z,pos_point_x,
     >     pos_point_y,pos_point_z,inside,particles_at_point,moat,
     >     number_particles,particles_real,pos_x,pos_y,pos_z,debug,
     >     tweaks)
      call timing(1,3)
c
      group=1
c
      call timing(-1,4)
      print*,'assign density'
      call assign_density(group,particles_at_point,
     >     hoc_particles,density,
     >     next_particles,grid_length,level_max,
     >     pos_x,pos_y,pos_z,particle_mass,pos_point_x,
     >     pos_point_y,pos_point_z,level,
     >     number_particles,periodic,hoc_points,next_points,
     >     particle_pointer,point_up_x,point_up_y,point_up_z,
     >     point_down_x,point_down_y,point_down_z,
     >     debug,point_pointer,inside,
     >     groups_maxx,points_maxx,particles_maxx,particles_real,
     >     tweaks)
c      print*,'assign density'
      call timing(1,4)
c     
      if(periodic) then
         print*,'periodic solver '
         if(debug)write(31,*)'periodic solver'
c     
         call timing(-1,5)
         call periodic_solver(grid_length,number_particles,
     >        density,potential_point,
     >        debug,groups_maxx,points_maxx,particles_maxx,tweaks)
         call timing(1,5)
      else
         if(debug)write(31,*)'isolated solver'
c     
         call timing(-1,6)
         call isolated_solver(grid_length,density,
     >        potential_point,
     >        debug,groups_maxx,points_maxx,particles_maxx,tweaks)
         call timing(1,6)
      end if
      if(debug)write(31,*)'force at point ',chain,group
c     
      call timing(-1,7)
      call force_at_point(group,hoc_points,next_points,
     >     force_point_x,force_point_y,force_point_z,
     >     potential_point,point_up_x,point_up_y,point_up_z,
     >     point_down_x,point_down_y,point_down_z,
     >     point_pointer,grid_length,level,inside,
     >     debug,groups_maxx,points_maxx,particles_maxx,tweaks)
      call timing(1,7)
c     
      work_to_be_done=.true.
c
      do while(work_to_be_done)
c     
         work_to_be_done=.false.
c     
         group_loop=hoc_groups
         do while(group_loop .gt. 0)
            group=group_loop
            call parts_in_group(group,hoc_points,next_points,
     >           hoc_particles,next_particles,particles_in_group)
c            print*,'what is going on ',group,level(group),level_max
            if(level(group) .eq. level_max) 
     >           number_high_groups(group)=0
            if(number_high_groups(group) .lt. 0) then
               work_to_be_done=.true.
c
               if(debug)write(31,*)'find high pairs',group,level(group)
c
               call timing(-1,8)
               call find_high_pairs(group,hoc_points,next_points,
     >              point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              particles_at_point,
     >              minimum_number,number_high_points,
     >              seq_to_list,high_pairs,
     >              list_high_1,list_high_2,
     >              debug,groups_maxx,points_maxx,
     >              points_maxx_group,particles_maxx,tweaks)
               call timing(1,8)
c     
               if(number_high_points .eq. 0) then
                  number_high_groups(group)=0
               else
                  if(debug)write(31,*)'equivalence class'
c
                  call timing(-1,9)
                  call equivalence_class(group_tmp,number_high_points,
     >                 list_high_1,list_high_2,high_pairs,
     >                 debug,groups_maxx,points_maxx,
     >                 points_maxx_group,particles_maxx,tweaks)
                  call timing(1,9)
c     
                  if(debug) then
                     do n=1,number_high_points
                        write(18,*)n,group_tmp(n)
                     end do
                  end if
                  mother_group_current=group
c     
                  if(debug)print*,' high groups'
c
                  call timing(-1,10)
                  call high_groups(highest_used_high_group,
     >                 number_high_groups,mother_group_current,
     >                 hoc_high_groups,next_high_groups,
     >                 number_high_points,hoc_high_points,
     >                 next_high_points,
     >                 mother_of_high_group,seq_to_list,group_tmp,
     >                 minimum_number,particles_at_point,
     >                 level,level_high,
     >                 point_up_x,point_up_y,point_up_z,point_down_x,
     >                 point_down_y,point_down_z,high_point_up_x,
     >                 high_point_up_y,high_point_up_z,
     >                 high_point_down_x,high_point_down_y,
     >                 high_point_down_z,
     >                 debug,groups_maxx,points_maxx,
     >                 high_points_maxx,particles_maxx,tweaks)
                  call timing(1,10)
c     
                  if(debug) print*,group,highest_used_group,
     >                 number_high_groups(group)
c
                  high_group=hoc_high_groups(mother_group_current)
                  do while(high_group .gt. 0)
                     new_group=highest_used_group+1
c
                     if(new_group .ge. groups_max) then
                        memory_value=3
                        print*,'memory value= ',memory_value
                        return
                     end if
c
                     next_groups(new_group)=hoc_groups
                     hoc_groups=new_group
                     mother_group(new_group)=group
                     if(debug)print*,' daughter group',group,
     >                    high_group,new_group,
     >                    number_high_groups(group),highest_used_point
c
                     call timing(-1,11)
                     call daughter_group(new_group,hoc_points,
     >                    next_points,hoc_high_points,next_high_points,
     >                    point_up_x,point_up_y,point_up_z,
     >                    point_down_x,point_down_y,point_down_z,
     >                    hoc_particles,next_particles,pos_x,pos_y,
     >                    pos_z,mother_group,high_group,inside,periodic,
     >                    level,pos_point_x,pos_point_y,pos_point_z,
     >                    point_pointer,level_max,high_point_up_x,
     >                    high_point_up_y,high_point_up_z,
     >                    particles_at_point,high_point_down_x,
     >                    high_point_down_y,high_point_down_z,
     >                    particle_pointer,
     >                    highest_used_point,highest_used_particle,
     >                    grid_length,debug,groups_maxx,points_maxx,
     >                    points_maxx_group,high_points_maxx,
     >                    particles_maxx,particles_real,memory_value,
     >                    tweaks)
                     if(memory_value .ne. 0) return
                     call timing(1,11)
c
                     number_high_groups(new_group)=-1
                     highest_used_group=new_group
c     
                     high_group=next_high_groups(high_group)
                  end do
               end if
            end if
            group_loop=next_groups(group_loop)
         end do
      end do
c     
      if(debug)write(31,*)'find highest level group'
c
      call timing(-1,12)
      call find_highest_level_group(
     >     hoc_groups,next_groups,
     >     hoc_points,next_points,hoc_particles,
     >     next_particles,particle_pointer,level,
     >     highest_level_group,debug,
     >     groups_maxx,points_maxx,particles_maxx,
     >     particles_real,tweaks)
      call timing(1,12)
c     
      group=1
      if(debug)write(31,*)'force at particles ',group
c     
      call timing(-1,13)
      call force_at_particle(group,
     >     hoc_points,next_points,potential_point,potential,
     >     force_point_x,force_point_y,force_point_z,
     >     particle_pointer,force_x,force_y,force_z,
     >     point_up_x,point_up_y,point_up_z,
     >     inside,hoc_particles,next_particles,
     >     pos_point_x,pos_point_y,
     >     pos_point_z,pos_x,pos_y,pos_z,level,
     >     level_max,grid_length,
     >     highest_level_group,debug,
     >     groups_maxx,points_maxx,
     >     particles_maxx,particles_real,tweaks)
      call timing(1,13)
c     
      if(debug)write(31,*)'find chains'
c
      number_chains=0
      call timing(-1,14)
      if(level_max .gt. 0) 
     >     call find_chains(hoc_groups,next_groups,mother_group,
     >     level,number_high_groups,number_chains,hoc_groups_chain,
     >     next_groups_chain,done_group,chain_started,chain_finished,
     >     groups_maxx,tweaks)
      call timing(1,14)
c
*** This finishes all the setup stuff, now we can get going.
c
      work_to_be_done=.true.
      do while(work_to_be_done)
         work_to_be_done=.false.
         do chain=1,number_chains
            group=hoc_groups_chain(chain)
            ok_start_chain(chain)=done_group(mother_group(group)) .and.
     >           .not. chain_finished(chain)
            do while(group .gt. 0 .and. ok_start_chain(chain))
               chain_started(chain)=.true.
               print*,'chain group ',chain,group,level(group)
               if(debug)write(31,*)'assign density ',chain,group
               call timing(-1,15)
               call assign_density(group,particles_at_point,
     >              hoc_particles,density,
     >              next_particles,grid_length,level_max,
     >              pos_x,pos_y,pos_z,particle_mass,pos_point_x,
     >              pos_point_y,pos_point_z,level,
     >              number_particles,periodic,hoc_points,next_points,
     >              particle_pointer,point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              debug,point_pointer,inside,
     >              groups_maxx,points_maxx,particles_maxx,
     >              particles_real,tweaks)
               call timing(1,15)
c     
               if(debug)write(31,*)'potential start ',chain,group
c     
               call timing(-1,16)
               call potential_start(group,hoc_points,next_points,
     >              potential_point,point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              point_pointer,grid_length,level,inside,
     >              debug,groups_maxx,points_maxx,particles_maxx,
     >              tweaks)
               call timing(1,16)
c     
               if(debug)write(31,*)'poisson solver ',chain,group
c     
               call timing(-1,17)
               print*,'poisson solver'
               call poisson_solver(group,hoc_points,next_points,
     >              point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              periodic,inside,potential_point,density,
     >              level,level_max,grid_length,
     >              debug,groups_maxx,points_maxx,particles_maxx,
     >              tweaks)
               call timing(1,17)
c     
               if(debug)write(31,*)'force at point ',chain,group
c     
               call timing(-1,18)
               call force_at_point(group,hoc_points,next_points,
     >              force_point_x,force_point_y,force_point_z,
     >              potential_point,point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              point_pointer,grid_length,level,inside,
     >              debug,groups_maxx,points_maxx,particles_maxx,
     >              tweaks)
               call timing(1,18)
c     
               if(debug)write(31,*)'force at particles ',chain,group
c     
               call timing(-1,19)
               call force_at_particle(group,
     >              hoc_points,next_points,potential_point,potential,
     >              force_point_x,force_point_y,force_point_z,
     >              particle_pointer,force_x,force_y,force_z,
     >              point_up_x,point_up_y,point_up_z,
     >              inside,hoc_particles,next_particles,
     >              pos_point_x,pos_point_y,
     >              pos_point_z,pos_x,pos_y,pos_z,level,
     >              level_max,grid_length,
     >              highest_level_group,debug,
     >              groups_maxx,points_maxx,
     >              particles_maxx,particles_real,tweaks)
               call timing(1,19)
c     
               print*,'finished chain group ',chain,group,level(group)
               done_group(group)=.true.
               group=next_groups_chain(group)
            end do
            print*,'finished chain ',chain
            chain_finished(chain)=chain_started(chain)
            work_to_be_done=.not. chain_finished(chain)
         end do
      end do
c     
      print*,'test it all'
      if(debug)write(31,*)'test it all'
      call timing(-1,20)
      if(debug .and. mod(tweaks/4,2) .eq. 1)
     >     call test_it_all(grid_length,level_max,hoc_groups,
     >     next_groups,hoc_points,next_points,inside,hoc_particles,
     >     next_particles,potential_point,potential,
     >     force_point_x,force_point_y,force_point_z,particle_pointer,
     >     pos_x,pos_y,pos_z,pos_point_x,pos_point_y,pos_point_z,
     >     point_pointer,force_x,force_y,force_z,
     >     level,highest_level_group,periodic,
     >     groups_maxx,points_maxx,particles_maxx,particles_real,
     >     tweaks)
      call timing(1,20)
c
      if(debug)write(31,*)'offset it back'
      call timing(-1,1)
      call offset(-1,random_offset,number_particles,
     >     grid_multiply,pos_x,pos_y,pos_z,tweaks)
      call timing(1,1)
c
      call timing(1,30)
      call timing(0,0)
c
      if(mod(tweaks/8,2) .eq. 1) then
         call dump_points(hoc_groups,next_groups,level,
     >        level_max,hoc_points,next_points,
     >        pos_point_x,pos_point_y,pos_point_z,
     >        pos_x,pos_y,pos_z,particle_pointer,hoc_particles,
     >        next_particles,
     >        groups_maxx,points_maxx,particles_maxx)
      end if
c     
      return
      end
