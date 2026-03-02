
      subroutine yfromfields(ny, u, v, f, y)

      implicit none
      integer ny
      double precision u(ny), v(ny), f(ny), y(ny)

      integer nymax
      include '../nymax.inc'
      double precision caux(2*nymax)
      integer i

      if (ny .gt. nymax) stop 'nymax too small in yfromfields'

      do i=1,ny
         caux(2*i-1) = u(i) / dble(ny)
         caux(2*i) = (v(i) + f(i)) / dble(ny)
      end do
      call four1(caux, ny, 1)
      call halve(caux, y, ny/2)

      return
      end
