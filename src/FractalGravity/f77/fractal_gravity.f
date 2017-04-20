      subroutine fractal_gravity(number_particles,grid_length,periodic,
     >     minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >     potential,force_x,force_y,force_z,
     >     highest_level_group,mother_group,
     >     padding,number_masks,level_mask,pos_mask,
     >     moat_0,random_offset,debug,tweaks,memory_value,
     >     parameters_integer,parameters_real,parameters_logical)

c     
c***  INPUT IS UNCHANGED EXCEPT PERIODIC DATA IS WRAPPED AROUND
c***  POSITIONS MAY BE OFFSET BUT THEY ARE PUT BACK BEFORE RETURNING
c     
      implicit none
c     
      include 'maxx.inc'
      include 'maxx_iso.inc'
      include 'maxes.inc'
c     
c     
c**** Memory_requirements=
c***  (L) 2  x  points_max
c***  (I) 14 x  points_max
c***  (R) 5  x  points_max
c***  (I) 2  x  particles_max
c***  (R) 9  x  number_particles
c***  (I) 1  x  number_particles
c***  (I) 2  x  maxx
c***  (R) 1 x grid_length_max**3
c***  plus assorted other small stuff
c     
c***  Input:
c***  number_particles (I)(Got that?)
c***  grid_length      (I)(number of cells in each dimension)
c***  periodic         (L) T/F is periodic/isolated boundary conditions
c***  minimum_number   (I) minimum number of particles in a cell so that 
c***  it is a high density cell
c***  level_max        (I) maximum level of groups
c***  pos_x            (R) x-positions for particles [0,1]
c***  pos_y            (R) y-positions for particles [0,1]
c***  pos_z            (R) z-positions for particles [0,1]
c***  if periodic, particles will be wrapped
c***  if isolated, outside particles will be ignored
c***  output:
c***  potential        (R) potential at particle
c***  force_x          (R) x-force at particles
c***  force_y          (R) y-force at particles
c***  force_z          (R) z-force at particles
c***  highest_level_group  (I) output only, ignored on input
c***  mother_group (I) output only, ignored on input
c***  padding (I)   pad each high_point with "padding" points in each direction
c***  number_of_mask (I) number of masks where we limit the levels
c***  level_mask (I) highest level for this mask
c***  pos_mask (R) [0,1] xmin,xmax,ymin,ymax,zmin,zmax for mask
c***  moat_0           (I) size of moat for isolated boundary conditions
c***  moat=min(moat_0,1). Not used in periodic calc
c***  random_offset    (I) if =0 no particle offset
c***  if != 0 seed for random offset
c***  debug            (L) T:debug run, F:no debug
c***  memory_value(I)  returned by daughter_group or other routines as error value
c***  if =0 no memory problem
c***  if =1 need more total points
c***  if =2 need more particles
c***  if =3 need more groups
c***  if =5 need more high points
c***  if =7 maxx too small
c***  if = 8 SOR did not converge
c***  anything else die
c***  tweaks
c***  mod(tweaks,2) .eq. 1 ; shrink masks with level of point
c***  mod(tweaks/2,2) .eq. 1 ; make buffer points
c***  mod(tweaks/4,2) .eq. 1 ; dump everything
c***  mod(tweaks/8,2) .eq. 1 ; use delta density
c***  mod(tweaks/16,2) .eq. 1 ; call make_rho_hat
c***  mod(tweaks/32,2) .eq. 1 ; call power_spectrum
c***  mod(tweaks/128,2) .eq. 1 ; use double precision SOR solver
c***  mod(tweaks/256,2) .eq. 1 ; smooth interface for forces 
c***  mod(tweaks/512,2) .eq. 1 ; smooth interface for density
c***  parameters_integer(1024) whatever I might later want to pass through to subroutines
c***  parameters_real(1024) whatever I might later want to pass through to subroutines
c***  parameters_logical(1024) whatever I might later want to pass through to subroutines
c***  parameters_integer(1) maxits_sor (maximum number of iterations for SOR)
c***  parameters_integer(2) spectrum_number (=0 powerlaw, =1 CDM)
c***  parameters_integer(3) level_max_zel (highest level generating power)
c***  parameters_integer(4) level_min (min level for splitting particles)
c***  parameters_integer(100+lev)(lev=0,level_max_zel) Initial seed for generating
c***  FFT of density filed at level lev.
c***  parameters_real(1) power_slope (slope of power law power spectrum or initial power law slope
c***  (if spectrum_number != 0)
c***  parameters_real(2) cutoff  (wave number of Gaussian cutoff of power spectrum)
c***  parameters_real(3) wave_0 (wave number for normalization of power spectrum)
c***  parameters_real(4) scaling (scaling in CDM spectrum)
c***  parameters_real(5) force_max (maximum force from one particle)
c***  parameters_real(6) ratio_rms_force (ratio_rms_force (transparent to user))
c***  parameters_real(7) ratio_rms_rho (ratio_rms_rho (transparent to user))
c***  parameters_real(8) epsilon_sor (accuracy parameter for SOR)
c***  parameters_logical(1) gauss_cut (is the the cutoff at wave_0 a Gaussian cutoff?) 
c***  green_0 is the density field with doubling in each dimension. In an attempt 
c***  to be clever and not waste space this array is stored in other arrays using equivalence statements.
c     
      integer number_particles,number_particles_tmp
      logical inside(points_max),periodic,work_to_be_done,debug
      logical done_group(groups_max),chain_finished(groups_max)
      logical chain_started(groups_max),ok_start_chain(groups_max)
      logical it_is_high(-1:points_max),heavy_version,delta_version
      logical parameters_logical(1024),start_up,force_smooth,pad_smooth
      logical small_group_allowed,buffer_it,dump_everything,shrink_mask
      logical double_solver,power_gen,set_zero,set_dens,set_scaling
c     
      real parameters_real(1024)
      real green(green_length),green_0(iso_maxx),epsilon_sor
      real pos_x(number_particles),pos_y(number_particles)
      real pos_z(number_particles),particle_mass(number_particles)
      real potential_point(points_max+1),density_0,sum_up
      real potential(number_particles),density(points_max)
      real force_point_x(points_max+1),force_point_y(points_max+1)
      real force_point_z(points_max)
      real force_x(number_particles),force_y(number_particles)
      real force_z(number_particles),force_max
      real pos_mask(20000)
      real sum_0(0:8),sum_1_x(0:8),sum_1_y(0:8),sum_1_z(0:8)
      real sum_2_x(0:8),sum_2_y(0:8),sum_2_z(0:8)
      real pos_x_tmp(particles_heavy),pos_y_tmp(particles_heavy)
      real pos_z_tmp(particles_heavy),particle_mass_tmp(particles_heavy)
      real edge_weight(points_max_edge)
c     
      integer parameters_integer(1024),level_min,maxits_sor
      integer random_offset,number_chains,memory_value,tweaks
      integer grid_length,level_max,level_max_zel,group,group_loop
      integer next_points(points_max),hoc_points(groups_max)
      integer hoc_particles(points_max),next_particles(particles_max)
      integer pos_point_x(points_max),pos_point_y(points_max)
      integer pos_point_z(points_max)
      integer point_up_x(points_max),point_up_y(points_max)
      integer point_up_z(points_max)
      integer point_down_x(points_max),point_down_y(points_max)
      integer point_down_z(points_max),point_pointer(points_max)
      integer particles_at_point(points_max)
      integer mother_of_high_group(groups_max)
      integer highest_used_point,highest_used_high_point
      integer highest_used_particle,highest_used_group,new_group
      integer particle_pointer(particles_max),real_pointer(points_max)
      integer lowest_point_in_group(groups_max)
      integer highest_point_in_group(groups_max)
      integer minimum_number,seq_to_list(high_points_max_group)
      integer high_pairs,point,points_maxx_edge
      integer list_high_1(3*high_points_max_group)
      integer list_high_2(3*high_points_max_group)
      integer group_tmp(high_points_max_group)
      integer next_groups(groups_max),hoc_groups
      integer level(groups_max),level_high(groups_max)
      integer number_high_groups(groups_max)
      integer highest_used_high_group,number_high_points
      integer mother_group_current,hoc_high_groups(groups_max)
      integer next_high_groups(groups_max),mother_group(groups_max)
      integer hoc_high_points(groups_max),particles_in_group(groups_max)
      integer points_in_group(groups_max)
      integer next_high_points(high_points_max)
      integer hoc_groups_chain(groups_max),next_groups_chain(groups_max)
      integer highest_level_group(number_particles),grid_multiply
      integer groups_maxx,points_maxx,particles_maxx,particles_real
      integer high_points_maxx,points_maxx_group,moat,moat_0,n
      integer high_group,chain,n_highs,particles_real_tmp
      integer really_high(high_points_max_group)
      integer padding,steps,p,lev,highest_level_init
      integer number_masks,level_mask(5000)
      integer particle_pointer_tmp(particles_max_0)
      integer hoc_particles_tmp(points_max),number_fakes
      integer next_particles_tmp(particles_heavy),particles_heavy_tmp
      integer highest_level_point(particles_max_0)
      integer fake_particle(points_max_edge)
c     
      equivalence (green_0(1),potential_point(1))
      equivalence (potential_point(points_max+1),force_point_x(1))
      equivalence (force_point_x(points_max+1),force_point_y(1))
      equivalence (force_point_y(points_max+1),force_point_z(1))
c     
      save green
      data steps/0/
c     
c***  groups_max(I) maximum number of groups
c***  points_max(I) maximum number of points
c***  high_points_max(I) maximum number of high points
c***  currently points_max=high_points_max
c***  points_max_group(I) maximum number of points in a single group
c***  particles_max(I) maximum number of particles
c     
c***  random_offset(I) random number seed
c***  if .ne. 0 offset particles to avoid interference with lattice
c***  if = 0 do nothing 
c***  inside(point;L) is this an inside point?
c***  periodic(L)  periodic boundary condition or isolate boundary conditions
c***  work_to_be_done (L) are there still groups possibly not done yet
c***  pos_x(particle;R) x position of particles 0< pos_x < 1
c***  pos_y(particle;R) y position of particles 0< pos_y < 1
c***  pos_z(particle;R) z position of particles 0< pos_z < 1
c***  potential_point(point;R) potential at a point
c***  potential(particle;R) potential at a particle
c***  density(point;R) density at a point
c***  force_point_x(point;R) x-force at a point
c***  force_point_y(point;R) y-force at a point
c***  force_point_z(point;R) z-force at a point
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
c***  0 <= pos_point_x < grid_length * 2**level_max
c***  pos_point_y(point;I) y-position of point 
c***  0 <= pos_point_y < grid_length * 2**level_max
c***  pos_point_z(point;I) z-position of point 
c***  0 <= pos_point_z < grid_length * 2**level_max
c***  point_up_x(point;I) neighbor point up in x (if < 0 no neigh point)
c***  point_up_y(point;I) neighbor point up in y (if < 0 no neigh point)
c***  point_up_z(point;I) neighbor point up in z (if < 0 no neigh point)
c***  point_down_x(point;I) neighbor point down in x (if < 0 no neigh point)
c***  point_down_y(point;I) neighbor point down in y (if < 0 no neigh point)
c***  point_down_z(point;I) neighbor point down in z (if < 0 no neigh point)
c***  particles_a_point(point;I) number of particles at point
c***  mother_of_high_group(high_group;I) mother group of high group
c***  highest_used_point(I) highest_used_point+1 is first new point available
c***  highest_used_high_point(I) 
c***  highest_used_high_point+1 is first new point available
c***  highest_used_particle(I)  
c***  highest_used_particle+1 first new virtual particle available
c***  highest_used_group(I)  
c***  highest_used_group+1 first new group available
c***  new group(I)   new daughter group just generated
c***  particle_pointer(particle;I) points to original particle in group 1
c***  minimum_number(I) minimum number of particles to make a high point 
c***  seq_to_list(high_point,I) pointer from seq list 
c***  of high points to link list of points
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
c***  of high points in group
c***  highest_level_group(particle,I) highest level group the 
c***  particle belongs to
c***  groups_maxx(I)    =groups_max   used for passing arrays
c***  points_maxx(I)    =points_max   used for passing arrays
c***  points_maxx_group(I)    =points_max_group   used for passing arrays
c***  particles_maxx(I) =particles_max   used for passing arrays
c***  particles_real(I) =number_particles   used for passing arrays
c***  moat_0(I) input buffer in isolated calculation
c***  moat(I)=min(moat_0,1)  buffer in isolated calculation
c***  n(I)        dummy variable
c     
c     
c     
c     
c***  Mr. Garibaldi to hapless Engineer: 
c***  It is a prototype, it does not have to be perfect,
c***  It just has to work.
c     
c***  Sheridan to Mr. Garibaldi:
c***  I am not interested in your problems, 
c***  I am only interested in your solutions.
c     
c***  Sirian Cybernetics Corporation
c***  Its fundamental design flaws are obscured
c***  by its superficial design flaws
c     
c***  This is a prototype of fractal gravity, it is not supposed to be pretty
c***  or fast, it is just supposed to work.
c***  This is a single node version
c     
c***  It "should" be simple to transform into a multi-node code.
c***  Code is called repeatedly, it does not remember anything
c***  from previous calls except a few useful arrays and the FT of the Green's
c***  function for an isolated simulation
c     
c***  JV Villumsen
c***  Should died many years ago.
c     
c***  call OFFSET
c***  offset postions by a small random shift, useful to avoid ringing between
c***  a uniform particle distribution and the uniform grid
c***  call NEGATIVES_GROUP
c***  set negative everything related to groups
c***  call NEGATIVES_POINTS
c***  set negative everything related to points
c***  call NEGATIVES_PARTICLES
c***  set negative everything related to particles
c***  call TREESTART
c***  generate group 1 at level 0.
c***  call ASSIGN_DENSITY
c***  find density field for group 1
c***  call PERIODIC_SOLVER
c***  call ISOLATED_SOLVER
c***  call periodic or isolated FFT solver for potential for group 1 points
c***  call FORCE_AT_POINT
c***  difference potential field to get force field for group 1 at points
c***  This finishes the group 1 force calculation.
c***  Now comes the tree generation
c***  For each group, find the high points, and arrange them into high groups.
c***  Each high group generates one daughter group. Continue this process
c***  recursively until we reach level=level_max or 
c***  there are no more high points.
c***  call HIGH_POINTS
c***  find all the high points in a group
c***  call BUFFER_POINTS
c***  if necessary, and if we so choose, add buffer high points
c***  call FIND_HIGH_PAIRS
c***  find pairs of high points
c***  call EQUIVALENCE_CLASS
c***  find equivalence classes of high points
c***  call HIGH_GROUPS
c***  make high groups from equivalence classes
c***  call DAUGHTER_GROUP
c***  for each high group generate a daughter group
c***  this finishes the tree of groups
c***  call FIND_HIGHEST_LEVEL_GROUP
c***  for each particle in the simulation, 
c***  find the groups it belongs to and choose the 
c***  one at the highest level
c***  call FORCE_AT_PARTICLE
C***  for all particles that exist only in group 1, 
c***  find the force on the particles
c***  call FIND_CHAINS
c***  generate chains of group useful for the potential solutions
c***  now loop over all the chains and the groups in a chain
c***  call ASSIGN DENSITY
c***  find the density at each point in the group
c***  call POTENTIAL_START
c***  use the potential values from the mother group to assign initial
c***  values to the potential at the group points
c***  call POISSON_SOLVER
c***  solve for the potential in the group 
c***  call FORCE_AT_POINT
c***  difference the potential to get the force field at points
c***  call FORCE_AT_PARTICLE
c***  interpolate to get the forces at the particles that belong to the group
c***  Now all is done
c***  call WRITE_IT
c***  print out everything
c***  if we so choose
c***  call OFFSET
c***  offset particles back
c***  THAT IS ALL FOR NOW FOLKS
c     
      call timing(-1,30)
c     
      steps=steps+1
      print*,'steps ',steps,tweaks
      if(number_particles .gt. particles_max_0) then
         print*,'number_particles ',number_particles,particles_max_0
         stop 'particles_max_0 not big enough'
      end if
c     
      if(periodic) then
         density_0=sum_up(particle_mass,1,number_particles,1)
      else
         density_0=0
      end if
c     
      small_group_allowed=.false.
      shrink_mask=mod(tweaks,2) .eq. 1
      buffer_it=mod(tweaks/2,2) .eq. 1
      dump_everything=mod(tweaks/4,2) .eq. 1
      delta_version=mod(tweaks/8,2) .eq. 1
      start_up=mod(tweaks/16,2) .eq. 1
      power_gen=mod(tweaks/32,2) .eq. 1
      double_solver=mod(tweaks/128,2) .eq. 1
      force_smooth=mod(tweaks/256,2) .eq. 1
      pad_smooth=mod(tweaks/512,2) .eq. 1
      number_fakes=0
      memory_value=0
      grid_multiply=grid_length*2**level_max
      groups_maxx=groups_max
      points_maxx=points_max
      high_points_maxx=high_points_max
      points_maxx_group=points_max_group
      points_maxx_edge=points_max_edge
      particles_maxx=particles_max
      particles_real=number_particles
      particles_heavy_tmp=particles_heavy
      maxits_sor=parameters_integer(1)
      level_min=parameters_integer(4)
      force_max=parameters_real(5)
      epsilon_sor=parameters_real(8)
      heavy_version=force_max .gt. 0 .and. delta_version
      if(.not. heavy_version) level_min=level_max+1
      if(.not. shrink_mask .and. (buffer_it .or. padding .gt. 0)) stop 
     >     'inconsistent masking'
c     
      it_is_high(-1)=.false.
      it_is_high(0)=.false.
c     
      moat=min(moat_0,1)
c     
      lowest_point_in_group(1)=0
c     
      if(g_maxx .lt. groups_maxx) then
         print*,'g_maxx too small ',g_maxx,groups_maxx
         memory_value=3
         return
      end if
c     
      if(periodic) then
         if((grid_length+1)**3 .gt. points_max) then
            memory_value=1
            print*,'case a memory value= ',memory_value,points_max,
     >           (grid_length+1)**3
            return
         end if
      end if
c     
      if(maxx .lt. points_max) then
         memory_value=7
         print*,'maxx too small'
         return
      end if
c     
      if(number_particles .gt. particles_max) then
         memory_value=2
         print*,'memory value= ',memory_value
         return
      end if
c     
      if(debug) write(41,*)'maxxes= ',
     >     particles_maxx,points_maxx,groups_maxx
c     
      if(debug)write(41,*)'offset it'
      call timing(-1,1)
      call offset(1,random_offset,number_particles,
     >     grid_multiply,pos_x,pos_y,pos_z,tweaks)
      call timing(1,1)
c     
      call timing(-1,2)

c     if(debug) then
c      write(41,*)'be negative'
      call negatives_groups(groups_maxx,hoc_points,
     >     mother_of_high_group,
     >     number_high_groups,next_groups,level,level_high,
     >     particles_in_group,points_in_group,
     >     hoc_high_groups,next_high_groups,mother_group,tweaks)
c     
      call negatives_points(points_maxx,high_points_maxx,
     >     points_maxx_group,next_points,hoc_particles,pos_point_x,
     >     pos_point_y,pos_point_z,point_up_x,point_down_x,point_up_y,
     >     point_down_y,point_up_z,point_down_z,inside,
     >     particles_at_point,density,
     >     potential_point,force_point_x,force_point_y,force_point_z,
     >     seq_to_list,
     >     list_high_1,list_high_2,group_tmp,next_high_points,tweaks)
c     end if
c     
      call negatives_particles(particles_maxx,particles_real,
     >     next_particles,
     >     particle_pointer,highest_level_group,
     >     potential,force_x,force_y,force_z,tweaks)
c     
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
      write(93,*)'step ',steps
      write(94,*)'step ',steps
c     
      call timing(-1,3)
      call tree_start(grid_length,level_max,group,groups_maxx,
     >     points_maxx,particles_maxx,periodic,hoc_points,
     >     next_points,highest_used_point,highest_used_particle,
     >     point_up_x,point_up_y,point_up_z,point_pointer,real_pointer,
     >     hoc_particles,next_particles,particle_pointer,
     >     point_down_x,point_down_y,point_down_z,pos_point_x,
     >     pos_point_y,pos_point_z,inside,particles_at_point,moat,
     >     number_particles,particles_real,pos_x,pos_y,pos_z,debug,
     >     tweaks)
      call timing(1,3)
c     
c******
      highest_point_in_group(1)=highest_used_point
c******
c     
      work_to_be_done=level_max .gt. 0
c     
      do while(work_to_be_done)
c     
         work_to_be_done=.false.
c     
         group_loop=hoc_groups
         do while(group_loop .gt. 0)
            group=group_loop
            if(debug)print*,'what is going on ',group,
     >           level(group),level_max
            if(level(group) .eq. level_max) 
     >           number_high_groups(group)=0
            if(number_high_groups(group) .lt. 0) then
               work_to_be_done=.true.
c     
               if(debug)write(41,*)'find high pairs',group,level(group)
c     
               call timing(-1,8)
c     
               if(debug)print*,'high points '
               call high_points(group,hoc_points,next_points,
     >              particles_at_point,minimum_number,
     >              really_high,n_highs,
     >              pos_point_x,pos_point_y,pos_point_z,number_masks,
     >              level_mask,pos_mask,grid_length,level_max,
     >              level,periodic,shrink_mask,
     >              it_is_high,groups_maxx,points_maxx,high_points_maxx,
     >              memory_value)
c     
               if(memory_value .ne. 0) return
c     
c     print*,'n_highs ',group,level(group),n_highs
               if(debug)print*,'buffer points '
               if(buffer_it .or. padding .gt.0) then
                  call buffer_points(group,really_high,n_highs,
     >                 it_is_high,level,grid_length,level_max,
     >                 pos_x,pos_y,pos_z,hoc_points,next_points,
     >                 pos_point_x,pos_point_y,pos_point_z,
     >                 point_up_x,point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 hoc_particles,
     >                 next_particles,particle_pointer,minimum_number,
     >                 periodic,inside,padding,
     >                 groups_maxx,points_maxx,
     >                 particles_maxx,particles_real)
               end if
c     
               call find_high_pairs(group,hoc_points,next_points,
     >              point_up_x,point_up_y,point_up_z,
     >              point_down_x,point_down_y,point_down_z,
     >              it_is_high,number_high_points,
     >              seq_to_list,high_pairs,
     >              list_high_1,list_high_2,
     >              debug,groups_maxx,points_maxx,
     >              points_maxx_group,particles_maxx,tweaks)
               call timing(1,8)
c     
               if(number_high_points .eq. 0) then
                  number_high_groups(group)=0
               else
                  if(debug)write(41,*)'equivalence class'
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
                  if(debug) print*,' high groups'
c     
                  call timing(-1,10)
                  call high_groups(highest_used_high_group,
     >                 number_high_groups,mother_group_current,
     >                 hoc_high_groups,next_high_groups,
     >                 number_high_points,hoc_high_points,
     >                 next_high_points,
     >                 mother_of_high_group,seq_to_list,group_tmp,
     >                 it_is_high,level,level_high,
     >                 point_up_x,point_up_y,point_up_z,point_down_x,
     >                 point_down_y,point_down_z,
     >                 debug,groups_maxx,points_maxx,
     >                 high_points_maxx,tweaks)
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
     >                    high_group,new_group,level(group),
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
     >                    point_pointer,level_max,it_is_high,
     >                    particles_at_point,particles_in_group,
     >                    points_in_group,particle_pointer,real_pointer,
     >                    lowest_point_in_group,highest_point_in_group,
     >                    highest_used_point,highest_used_particle,
     >                    grid_length,debug,groups_maxx,points_maxx,
     >                    points_maxx_group,high_points_maxx,
     >                    particles_maxx,particles_real,
     >                    memory_value,tweaks)
                     if(memory_value .ne. 0) return
                     call timing(1,11)
c     
                     if(points_in_group(new_group) .lt. 27) stop 
     >                    'bad new group'
                     if(points_in_group(new_group) .eq. 27 .and.
     >                    .not. small_group_allowed) then
                        number_high_groups(new_group)=0
c***  This means do not generate any groups from new_group
                     else
                        number_high_groups(new_group)=-1
                     end if
c     
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
      if(debug)write(41,*)'find highest level group'
c     
      call timing(-1,12)
      call find_highest_level_group(
     >     hoc_groups,next_groups,
     >     hoc_points,next_points,hoc_particles,
     >     next_particles,particle_pointer,number_particles,
     >     heavy_version,force_max,mother_group,
     >     particle_mass,grid_length,
     >     point_down_x,point_down_y,point_down_z,
     >     real_pointer,point_pointer,highest_level_point,
     >     level,highest_level_group,debug,
     >     groups_maxx,points_maxx,particles_maxx,
     >     particles_real,tweaks)
      call timing(1,12)
c     
      if(force_smooth) then
         call where_is_edge(highest_level_group,mother_group,
     >        highest_level_point,number_particles,
     >        fake_particle,point_pointer,real_pointer,
     >        point_up_x,point_up_y,point_up_z,
     >        point_down_x,point_down_y,point_down_z,inside,
     >        number_fakes,groups_maxx,
     >        points_maxx,points_maxx_edge,
     >        particles_max_0)
      end if
c     

c     do p=1,number_particles
c     if(mod(p,100) .eq. 1) print*,'pp ',p,highest_level_group(p)
c     end do
c     
      if(debug)write(41,*)'find chains'
c     
      number_chains=0
      call timing(-1,14)
      if(level_max .gt. 0) 
     >     call find_chains(hoc_groups,next_groups,mother_group,
     >     level,number_high_groups,number_chains,hoc_groups_chain,
     >     next_groups_chain,done_group,chain_started,chain_finished,
     >     debug,groups_maxx,tweaks)
      call timing(1,14)
c     
***   This finishes all the setup stuff, now we can get going.
      group=1
c     
      call timing(-1,4)
      if(debug)print*,'assign density'
      set_zero=.true.
      set_dens=.true.
      set_scaling=.true.
      call assign_density(group,hoc_particles,density,
     >     next_particles,grid_length,level_max,
     >     pos_x,pos_y,pos_z,particle_mass,density_0,pos_point_x,
     >     pos_point_y,pos_point_z,level,
     >     number_particles,periodic,hoc_points,next_points,
     >     particle_pointer,point_up_x,point_up_y,point_up_z,
     >     point_down_x,point_down_y,point_down_z,
     >     debug,point_pointer,inside,set_zero,set_dens,set_scaling,
     >     groups_maxx,points_maxx,particles_maxx,particles_real,
     >     tweaks)
      if(debug)print*,'assign density'
      call timing(1,4)
c     
      if(periodic) then
         if(debug)print*,'periodic solver '
         if(debug)write(41,*)'periodic solver'
c     
         call timing(-1,5)
         call periodic_solver(grid_length,number_particles,
     >        density,potential_point,0,
     >        debug,groups_maxx,points_maxx,particles_maxx,tweaks,
     >        parameters_integer,parameters_real,parameters_logical)
         call timing(1,5)
      else
         if(debug)write(41,*)'isolated solver'
c     
         call timing(-1,6)
         write(63,*)'iso a ',grid_length,debug,memory_value,
     >        groups_maxx,points_maxx,particles_maxx,tweaks
c     
         if(iso_maxx .gt. 4*points_max) then
            memory_value=2
            print*,'isolated length ',iso_maxx,4*points_max
            return
         end if
c     
         call isolated_solver(grid_length,density,
     >        potential_point,debug,memory_value,
     >        green,green_0,
     >        groups_maxx,points_maxx,particles_maxx,tweaks)
         write(63,*)'iso b ',grid_length,debug,memory_value,
     >        groups_maxx,points_maxx,particles_maxx,tweaks
c     
         call timing(1,6)
         if(memory_value .ne. 0) return
      end if
      if(debug)write(41,*)'force at point ',chain,group
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
      if(debug)write(41,*)'force at particles ',group
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
      level_max_zel=parameters_integer(3)
      highest_level_init=-1
c     
      if(level_max .gt. 0 .and. level_max_zel .gt. 0 .and. 
     >     start_up) then
         group=hoc_groups
         do while(group .gt. 0)
            highest_level_init=max(highest_level_init,level(group))
c     
            group=next_groups(group)
         end do
      end if 
c     
      highest_level_init=min(highest_level_init,level_max_zel,level_max)
c     
      lev=1
      do while(lev .le. level_max .and. start_up)
         if(highest_used_point+grid_length**3 .gt. points_maxx)then
            print*,'not enough points to generate initial conditions ',
     >           highest_used_point,grid_length**3,
     >           highest_used_point+grid_length**3
            stop 'I need more points'
         end if
         if(lev .le. highest_level_init) then
c     
            call periodic_solver(grid_length,number_particles,
     >           density,potential_point(highest_used_point+1),lev,
     >           debug,groups_maxx,points_maxx,particles_maxx,tweaks,
     >           parameters_integer,parameters_real,parameters_logical)
         end if
c     
         group=hoc_groups 
         do while(group .gt. 0)
            if(level(group) .eq. lev) then
               call assign_force_pot(potential_point,highest_used_point,
     >              hoc_points(group),next_points,
     >              force_point_x,force_point_y,force_point_z,
     >              pos_point_x,pos_point_y,pos_point_z,lev,
     >              highest_level_init,level_max,
     >              grid_length,points_maxx,tweaks,
     >              parameters_integer,parameters_real,
     >              parameters_logical)
            end if
            group=next_groups(group)
         end do
         lev=lev+1
      end do
c     
      work_to_be_done=number_chains .gt. 0
c     
      if(work_to_be_done) then
c     
         number_particles_tmp=0
c     
         if(heavy_version) then
            call heavies(heavy_version,grid_length,pos_x,pos_y,pos_z,
     >           pos_x_tmp,pos_y_tmp,pos_z_tmp,
     >           particle_mass,particle_mass_tmp,
     >           level,number_particles,number_particles_tmp,
     >           highest_level_group,level_min,level_max,force_max,
     >           particles_heavy_tmp,groups_maxx,particles_real)
c     
         end if
c     
         if((buffer_it .or. padding .gt. 0) .and. pad_smooth) then
            call edge_fakes(grid_length,pos_x,pos_y,pos_z,
     >           pos_x_tmp,pos_y_tmp,pos_z_tmp,
     >           particle_mass,particle_mass_tmp,
     >           level,number_particles,number_particles_tmp,
     >           highest_level_point,point_down_x,point_down_y,
     >           point_down_z,it_is_high,particles_at_point,
     >           point_pointer,real_pointer,minimum_number,
     >           highest_level_group,level_min,level_max,
     >           particles_heavy_tmp,groups_maxx,points_maxx,
     >           particles_real)
         end if
c     
         particles_real_tmp=number_particles_tmp
c     
         print*,' starting lists ',grid_length,level_max,hoc_groups,
     >        number_particles_tmp,groups_maxx,points_maxx,
     >        particles_maxx,particles_real
         if(number_particles_tmp .gt. 0)
     >        call particle_lists(grid_length,level_max,level,
     >        hoc_groups,next_groups,point_pointer,
     >        pos_x_tmp,pos_y_tmp,pos_z_tmp,
     >        particle_pointer_tmp,hoc_particles_tmp,next_particles_tmp,
     >        pos_point_x,pos_point_y,pos_point_z,
     >        point_up_x,point_up_y,point_up_z,
     >        hoc_points,next_points,real_pointer,
     >        periodic,number_particles_tmp,
     >        groups_maxx,points_maxx,particles_maxx,particles_real_tmp)
c     
         print*,' ending lists'
      end if
c     
      do while(work_to_be_done)
         work_to_be_done=.false.
         do chain=1,number_chains
            group=hoc_groups_chain(chain)
            ok_start_chain(chain)=done_group(mother_group(group)) .and.
     >           .not. chain_finished(chain)
            do while(group .gt. 0 .and. ok_start_chain(chain))
               chain_started(chain)=.true.
               if(debug)print*,'chain group ',chain,group,level(group)
               if(debug)write(41,*)'assign density ',chain,group
c**** 
               if(.not. start_up) then
                  call timing(-1,15)
c     
                  set_zero=.true.
                  set_dens=number_particles_tmp .le. 0
                  set_scaling= number_particles_tmp .le. 0
                  call assign_density(group,hoc_particles,density,
     >                 next_particles,grid_length,level_max,
     >                 pos_x,pos_y,pos_z,particle_mass,density_0,
     >                 pos_point_x,pos_point_y,pos_point_z,level,
     >                 number_particles,periodic,hoc_points,next_points,
     >                 particle_pointer,point_up_x,
     >                 point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 debug,point_pointer,inside,
     >                 set_zero,set_dens,set_scaling,
     >                 groups_maxx,points_maxx,particles_maxx,
     >                 particles_real,tweaks)
c     
                  if(number_particles_tmp .gt. 0) then
                     set_zero=.false.
                     set_dens=.true.
                     set_scaling=.true.
                     call assign_density(group,hoc_particles_tmp,
     >                    density,
     >                    next_particles_tmp,grid_length,level_max,
     >                    pos_x_tmp,pos_y_tmp,pos_z_tmp,
     >                    particle_mass_tmp,density_0,
     >                    pos_point_x,pos_point_y,pos_point_z,level,
     >                    number_particles_tmp,periodic,
     >                    hoc_points,next_points,
     >                    particle_pointer_tmp,point_up_x,
     >                    point_up_y,point_up_z,
     >                    point_down_x,point_down_y,point_down_z,
     >                    debug,point_pointer,inside,
     >                    set_zero,set_dens,set_scaling,
     >                    groups_maxx,points_maxx,particles_maxx,
     >                    particles_real_tmp,tweaks)
                  end if
                  call timing(1,15)
c     
                  if(debug)write(41,*)'potential start ',chain,group
c     
                  call timing(-1,16)
c     
                  call potential_start(group,hoc_points,next_points,
     >                 potential_point,point_up_x,point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 point_pointer,grid_length,level,inside,
     >                 debug,groups_maxx,points_maxx,particles_maxx,
     >                 tweaks)
                  call timing(1,16)
c     
                  if(debug)write(41,*)'poisson solver ',chain,group
c     
                  if(debug)print*,'poisson solver'
c     
                  call timing(-1,17)
                  call poisson_solver(group,hoc_points,next_points,
     >                 point_up_x,point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 periodic,inside,potential_point,density,
     >                 level,level_max,grid_length,
     >                 maxits_sor,epsilon_sor,memory_value,
     >                 debug,groups_maxx,points_maxx,particles_maxx,
     >                 double_solver)
c     
                  call timing(1,17)
                  if(memory_value .ne. 0) return
               end if
c     
               if(debug)write(41,*)'force at point ',chain,group
c     
               call timing(-1,18)
c     
               if(.not. start_up) then
                  call force_at_point(group,hoc_points,next_points,
     >                 force_point_x,force_point_y,force_point_z,
     >                 potential_point,point_up_x,point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 point_pointer,grid_length,level,inside,
     >                 debug,groups_maxx,points_maxx,particles_maxx,
     >                 tweaks)
               else
c     
                  call force_at_point_zel(group,hoc_points,next_points,
     >                 force_point_x,force_point_y,force_point_z,
     >                 point_up_x,point_up_y,point_up_z,
     >                 point_down_x,point_down_y,point_down_z,
     >                 point_pointer,grid_length,level,
     >                 highest_level_init,inside,
     >                 debug,groups_maxx,points_maxx,particles_maxx,
     >                 tweaks)
               end if
               call timing(1,18)
c     
               if(debug)write(41,*)'force at particles ',chain,group
c     
               if(debug)print*,'finished chain group ',
     >              chain,group,level(group)
               done_group(group)=.true.
               group=next_groups_chain(group)
            end do
            if(debug)print*,'finished chain ',chain
            chain_finished(chain)=chain_started(chain)
            work_to_be_done=.not. chain_finished(chain)
         end do
      end do
c     
      call timing(-1,19)
c     
      group=hoc_groups
      do while(group .gt. 0)
         if(group .gt. 1) then
c     
            call force_at_particle(group,
     >           hoc_points,next_points,
     >           potential_point,potential,
     >           force_point_x,force_point_y,force_point_z,
     >           particle_pointer,force_x,force_y,force_z,
     >           point_up_x,point_up_y,point_up_z,
     >           inside,hoc_particles,next_particles,
     >           pos_point_x,pos_point_y,
     >           pos_point_z,pos_x,pos_y,pos_z,level,
     >           level_max,grid_length,
     >           highest_level_group,debug,
     >           groups_maxx,points_maxx,
     >           particles_maxx,particles_real,tweaks)
            
         end if
         
         group=next_groups(group)
      end do
      call timing(1,19)
      
      if(force_smooth .and. number_fakes .gt. 0) then
         call timing(-1,19)
c     
         print*,'fakes ',steps,number_fakes
c     
         call edge_distance(fake_particle,
     >        highest_level_group,highest_level_point,
     >        number_fakes,pos_x,pos_y,pos_z,pos_point_x,pos_point_y,
     >        pos_point_z,point_up_x,point_up_y,point_up_z,
     >        point_down_x,point_down_y,point_down_z,
     >        periodic,inside,grid_length,level,level_max,real_pointer,
     >        edge_weight,groups_maxx,particles_real,
     >        points_maxx,points_maxx_group,particles_maxx,tweaks)
c     
         call force_at_fake_particles(potential_point,potential,
     >        force_point_x,force_point_y,force_point_z,
     >        force_x,force_y,force_z,
     >        point_up_x,point_up_y,point_up_z,
     >        point_down_x,point_down_y,point_down_z,
     >        pos_point_x,pos_point_y,pos_point_z,
     >        pos_x,pos_y,pos_z,level,level_max,grid_length,edge_weight,
     >        real_pointer,point_pointer,mother_group,
     >        highest_level_group,highest_level_point,
     >        fake_particle,number_fakes,debug,
     >        groups_maxx,points_maxx,particles_maxx,particles_real,
     >        tweaks)
         
         call timing(1,19)
      end if
c     
      call timing(-1,20)
      if(dump_everything) then
         write(97)highest_used_group,highest_used_point,
     >        highest_used_particle,number_particles,
     >        hoc_groups,
     >        groups_max,points_max,particles_max
         write(98)highest_used_group,highest_used_point,
     >        highest_used_particle,number_particles,
     >        hoc_groups,
     >        groups_max,points_max,particles_max
         write(99)highest_used_group,highest_used_point,
     >        highest_used_particle,number_particles,
     >        hoc_groups,
     >        groups_max,points_max,particles_max
         call write_int(97,next_groups,highest_used_group)
         call write_int(97,hoc_points,highest_used_group)
         call write_int(97,hoc_high_points,highest_used_group)
         call write_int(97,lowest_point_in_group,highest_used_group)
         call write_int(97,highest_point_in_group,highest_used_group)
         call write_int(97,level,highest_used_group)
         call write_int(97,number_high_groups,highest_used_group)
         call write_int(97,mother_group,highest_used_group)
         call write_int(97,particles_in_group,highest_used_group)
         call write_int(97,points_in_group,highest_used_group)
         call write_int(97,hoc_groups_chain,highest_used_group)
         call write_int(97,next_groups_chain,highest_used_group)
c     
         call write_int(98,next_points,highest_used_point)
         call write_int(98,hoc_particles,highest_used_point)
         call write_int(98,pos_point_x,highest_used_point)
         call write_int(98,pos_point_y,highest_used_point)
         call write_int(98,pos_point_z,highest_used_point)
         call write_int(98,point_up_x,highest_used_point)
         call write_int(98,point_up_y,highest_used_point)
         call write_int(98,point_up_z,highest_used_point)
         call write_int(98,point_down_x,highest_used_point)
         call write_int(98,point_down_y,highest_used_point)
         call write_int(98,point_down_z,highest_used_point)
         call write_int(98,point_pointer,highest_used_point)
         call write_int(98,particles_at_point,highest_used_point)
c     
         call write_log(98,inside,highest_used_point)
         call write_log(98,it_is_high(1),highest_used_point)
c     
         call write_real(98,density,highest_used_point)
c*******
         call write_real(98,density,highest_used_point)
c*******
         call write_real(98,potential_point,highest_used_point)
         call write_real(98,force_point_x,highest_used_point)
         call write_real(98,force_point_y,highest_used_point)
         call write_real(98,force_point_z,highest_used_point)
c     
         call write_int(99,next_particles,highest_used_particle)
         call write_int(99,particle_pointer,highest_used_particle)
c     
         call write_real(99,pos_x,number_particles)
         call write_real(99,pos_y,number_particles)
         call write_real(99,pos_z,number_particles)
         call write_real(99,potential,number_particles)
         call write_real(99,force_x,number_particles)
         call write_real(99,force_y,number_particles)
         call write_real(99,force_z,number_particles)
         call write_real(99,particle_mass,number_particles)
c     
         call write_int(99,highest_level_group,number_particles)

      end if
      call timing(1,20)
c     
      if(debug)write(41,*)'offset it back'
      call timing(-1,1)
      call offset(-1,random_offset,number_particles,
     >     grid_multiply,pos_x,pos_y,pos_z,tweaks)
      call timing(1,1)
c     
      call timing(1,30)
      call timing(0,0)
c     
      if((buffer_it .or. padding .gt. 0) .and. debug) then
         group=hoc_groups
         do while(group .gt. 0)
            point=hoc_points(group)
            do while(point .gt. 0)
               if(point_pointer(point) .gt. 0) then
                  if(.not. inside(point_pointer(point))) then
                     if(mother_group(group) .eq. 1 .and. 
     >                    .not. periodic) then
                        print*,'touching moat ',group,point,
     >                       point_pointer(point)
                     else
                        print*,'buffer not perfet ',group,
     >                       mother_group(group),level(group),
     >                       point,point_pointer(point),
     >                       pos_point_x(point),
     >                       pos_point_y(point),
     >                       pos_point_z(point)
                     end if
                  end if
               end if
               point=next_points(point)
            end do
            group=next_groups(group)
         end do
      end if
c     
      do lev=0,level_max
         sum_0(lev)=1.0e-10
         sum_1_x(lev)=0.0
         sum_1_y(lev)=0.0
         sum_1_z(lev)=0.0
         sum_2_x(lev)=0.0
         sum_2_y(lev)=0.0
         sum_2_z(lev)=0.0
      end do
c     
      do p=1,number_particles
         lev=level(highest_level_group(p))
         sum_0(lev)=sum_0(lev)+1.0
         sum_1_x(lev)=sum_1_x(lev)+force_x(p)
         sum_1_y(lev)=sum_1_y(lev)+force_y(p)
         sum_1_z(lev)=sum_1_z(lev)+force_z(p)
         sum_2_x(lev)=sum_2_x(lev)+force_x(p)**2
         sum_2_y(lev)=sum_2_y(lev)+force_y(p)**2
         sum_2_z(lev)=sum_2_z(lev)+force_z(p)**2
      end do
c     
      do lev=0,level_max
         write(68,68)steps,lev,sum_0(lev),
     >        sum_1_x(lev)/sum_0(lev),sum_1_y(lev)/sum_0(lev),
     >        sum_1_z(lev)/sum_0(lev),
     >        sum_2_x(lev)/sum_0(lev),sum_2_y(lev)/sum_0(lev),
     >        sum_2_z(lev)/sum_0(lev)
 68      format(i5,i3,7(1pe13.5))
      end do
c     
      write(94,97)highest_used_group,highest_used_particle,
     >     highest_used_point
 97   format(' highest group,particle,point ',3i8)
      return
      end
