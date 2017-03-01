      subroutine edge_distance(fake_particle,
     >     highest_level_group,highest_level_point,
     >     number_fakes,pos_x,pos_y,pos_z,pos_point_x,pos_point_y,
     >     pos_point_z,up_x,up_y,up_z,down_x,down_y,down_z,
     >     periodic,inside,grid_length,level,level_max,real_pointer,
     >     edge_weight,groupies,particles_real,long,long_edge,
     >     particles_maxx,tweaks)
c     
      implicit none
c     
      integer groupies,long,long_edge,particles_maxx,particles_real
      integer tweaks,number_fakes,rp,rp1,rp2,rp3
c     
      logical inside(long),edge_list,ok,periodic
c     
      real edge_weight(number_fakes)
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real)
      real scaling,x,y,z,dist,i,j,pabs
c     
      integer point,particle,p,edges,list(8)
      integer level(groupies)
      integer up_x(long),up_y(long)
      integer up_z(long)
      integer down_x(long),down_y(long)
      integer down_z(long)
      integer part,real_pointer(long)
      integer pos_point_x(long),pos_point_y(long)
      integer pos_point_z(long)
      integer grid_length,level_max,highest_level_group(particles_real)
      integer highest_level_point(particles_real)
      integer group,n,grid_multiply
      integer fake_particle(number_fakes)
c     
      pabs(i,j)=min(abs(i),abs(i-j),abs(i+j))
c
      grid_multiply=grid_length*2**level_max
      scaling=grid_multiply
c     
      do part=1,number_fakes
         particle=fake_particle(part)
         point=highest_level_point(particle)
         group=highest_level_group(particle)
         x=pos_x(particle)*scaling
         y=pos_y(particle)*scaling
         z=pos_z(particle)*scaling
c     
         p=point
         rp=real_pointer(point)-1
         rp3=rp/9
         rp2=mod(rp/3,3)
         rp1=mod(rp,3)
         if(rp3 .eq. 1) p=down_z(p)
         if(rp2 .eq. 1) p=down_y(p)
         if(rp1 .eq. 1) p=down_x(p)
         ok=edge_list(p,inside,list,edges,
     >        up_x,up_y,up_z,down_x,down_y,down_z,long)
c
         if(.not. ok) stop 'badd'
c     
c         divide=2.0**(level_max-level(group))
         dist=-1.0
         do n=1,edges
            if(periodic) then
               dist=max(dist,
     >              pabs(x-float(pos_point_x(list(n))),grid_multiply),
     >              pabs(y-float(pos_point_y(list(n))),grid_multiply),
     >              pabs(z-float(pos_point_z(list(n))),grid_multiply))
            else
               dist=max(dist,
     >              abs(x-float(pos_point_x(list(n)))),
     >              abs(y-float(pos_point_y(list(n)))),
     >              abs(z-float(pos_point_z(list(n)))))
            end if
         end do
         dist=dist*2.0**(level(group)-level_max)
         edge_weight(part)=dist-1.0
         if(abs(edge_weight(part)-0.5) .gt. 0.5) then
         write(95,95)part,particle,dist,edge_weight(part),
     >        group,point,level(group),2**(-level(group)+level_max)
 95      format(2i8,2(1pe13.5),4i8)
         stop 'overweight'
      end if
      end do
      return
      end
c
      logical function edge_list(p,ins,list,length,
     >     up_x,up_y,up_z,down_x,down_y,down_z,long)
c
      implicit none
c
      integer long
      logical ins(long)
      integer length,p,list(8)
      integer up_x(long),up_y(long),up_z(long)
      integer down_x(long),down_y(long),down_z(long)
c
      length=0
c
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=up_x(up_x(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=up_y(up_y(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=down_x(down_x(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=up_z(up_z(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=up_x(up_x(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=down_y(down_y(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
      p=down_x(down_x(p))
      if(.not. ins(p)) then
         length=length+1
         list(length)=p
      end if
c
      edge_list=length .gt. 0
      return
      end
