
      subroutine mprove(a, alud, n, np, indx, b, x)

      integer n, np, indx(n), nmax
      double precision a(np,np), alud(np,np), b(n), x(n)
      parameter (nmax=1024)
      integer i, j
      double precision r(nmax)
      double precision sdp, maximum

      do i=1,n
         sdp = -b(i)
         do j=1,n
            sdp = sdp + a(i,j) * x(j)
         end do
         r(i) = sdp
      end do
      call lubksb(alud, n, np, indx, r)
      maximum = 0.d0
      do i=1,n
         maximum = max( abs(r(i)), maximum )
         x(i) = x(i) - r(i)
      end do
      write(*,*) 'Max(mprove)=', maximum

      return

      end
