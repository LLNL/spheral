      real function sum_up(x,n_1,n_2,n_d)
c     
      implicit none
c     
      real x(*)
      integer n,n_1,n_2,n_d
c     
      sum_up=0.0
      do n=n_1,n_2,n_d
       sum_up=sum_up+x(n)
      end do
c     
      return
      end
