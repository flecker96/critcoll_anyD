      subroutine double(x, xx, n)
      implicit none
      integer n, i
C     x and xx can safely be the same array (of size 4*n or bigger)

      double precision x(2*n), xx(4*n)
      do i=1,n
         xx(i) = x(i)
      end do
      do i=n+1,2*n
         xx(2*n+i) = x(i)
      end do
      xx(n+1) = 0.5d0 * xx(3*n+1)
      xx(3*n+1) = 0.5d0 * xx(3*n+1)
      xx(n+2) = 0.5d0 * xx(3*n+2)
      xx(3*n+2) = 0.5d0 * xx(3*n+2)
      do i=n+3,3*n
         xx(i) = 0.d0
      end do      

      return
      end

