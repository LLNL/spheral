      real function sum_prod(n1,n2,n3,x,y)
c     
      implicit none
c     
      integer n,n1,n2,n3
      real x(*),y(*)
c     
      sum_prod=0.0
      do n=n1,n2,n3
       sum_prod=sum_prod+x(n)*y(n)
      end do
c     
      return
      end
