      subroutine set_constant_real(x,n_1,n_2,n_d,what)
c     
      implicit none
c     
      integer n,n_1,n_2,n_d
      real x(n_2),what
c     
      do n=n_1,n_2,n_d
       x(n)=what
      end do
c     

      return
      end
