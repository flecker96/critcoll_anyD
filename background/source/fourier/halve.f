      subroutine halve(xx, x, n)
      implicit none
      integer n, i
      double precision xx(4*n), x(2*n)
C     xx and x can safely be the same array (of size 4*n or bigger)

      do i=1,n
         x(i) = xx(i)
      end do
      x(n+1) = xx(n+1) + xx(3*n+1)
      x(n+2) = xx(n+2) + xx(3*n+2)
      do i=n+3,2*n
         x(i) = xx(2*n+i)
      end do

      return
      end
