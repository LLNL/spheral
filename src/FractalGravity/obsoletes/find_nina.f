      subroutine find_nina(point,adj,
     >  up_x,up_y,up_z,down_x,down_y,down_z,n_sym,points_maxx)
c     
      implicit none
c     
      integer points_maxx
      integer point,i,j,k,up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx)
      integer down_y(points_maxx),down_z(points_maxx)
      integer adj(-1:1,-1:1,-1:1),n,g,n_plus,n_plus_0
      integer adj_tmp(-4:4,-4:4,-4:4),n_sym
c     
      do k=-n_sym,n_sym
        do j=-n_sym,n_sym
          do i=-n_sym,n_sym
            adj_tmp(i,j,k)=-1
          end do
        end do
      end do
c     
      adj_tmp(0,0,0)=point
c     
      n_plus_0=-1
      n_plus=0
      do n=1,27
        if(n_plus .gt. n_plus_0) then
          n_plus_0=n_plus
          n_plus=0
          do k=-n_sym,n_sym
            do j=-n_sym,n_sym
              do i=-n_sym,n_sym
                g=adj_tmp(i,j,k)
                if(g .gt. 0) then
                  n_plus=n_plus+1
                  if(k .lt. n_sym) adj_tmp(i,j,k+1)=up_z(g)
                  if(k .gt. -n_sym) adj_tmp(i,j,k-1)=down_z(g)
                  if(j .lt. n_sym) adj_tmp(i,j+1,k)=up_y(g)
                  if(j .gt. -n_sym) adj_tmp(i,j-1,k)=down_y(g)
                  if(i .lt. n_sym) adj_tmp(i+1,j,k)=up_x(g)
                  if(i .gt. -n_sym) adj_tmp(i-1,j,k)=down_x(g)
                end if
              end do
            end do
          end do
        end if
      end do
c
      do k=-1,1
        do j=-1,1
          do i=-1,1
            adj(i,j,k)=adj_tmp(i,j,k)
          end do
        end do
      end do
c     
      return
      end
