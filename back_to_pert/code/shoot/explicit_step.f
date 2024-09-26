C     nu-stage Implicit Runge-Kutta (2nu-order)

      subroutine explicit_step(ny, xin, yin, xout, yout,
     $     derivs, nu, crit, maxits, itsreach,
     $     junkc, junkf1, junkf2)

      implicit none
      integer ny, nu, maxits, itsreach
      double precision xin, yin(ny), xout, yout(ny), crit
      external derivs
      double precision junkc, junkf1(ny), junkf2(ny)

      integer nymax, numax
      include '../nymax.inc'

      parameter(numax=4)
      double precision xi(nymax,numax), xiold(nymax,numax), sum
      double precision dx, x(numax), f(nymax,numax)
      double precision a(numax,numax), b(numax), c(numax)

      integer i, ii, j, its

      if (ny .gt. nymax) stop 'nymax too small in implicit_step'

C     Define some ERK nu=2,3,4. Note that they are not unique.
C     Initialize everything to zero.
      do i=1,nu
         do ii=1,nu
            a(i,ii) = 0.d0
         end do
         b(i) = 0.d0
         c(i) = 0.d0
      end do
      if (nu.eq.1) then
         write(*,*) 'The simple Euler rule is unstable'
         stop
      else if (nu.eq.2) then
         a(2,1) = 0.5d0
         b(2) = 1.d0
         c(2) = 0.5d0
      else if (nu.eq.3) then
         a(2,1) = 0.5d0
         a(3,1) = -1.d0
         a(3,2) = 2.d0
         b(1) = 1.d0/6.d0
         b(2) = 2.d0/3.d0
         b(3) = 1.d0/6.d0
         c(2) = 0.5d0
         c(3) = 1.d0
      else if (nu.eq.4) then
         a(2,1) = 0.5d0
         a(3,2) = 0.5d0
         a(4,3) = 1.d0
         b(1) = 1.d0/6.d0
         b(2) = 1.d0/3.d0
         b(3) = 1.d0/3.d0
         b(4) = 1.d0/6.d0
         c(2) = 0.5d0
         c(3) = 0.5d0
         c(4) = 1.d0
      else
         write(*,*) 'ERK nu=', nu, ' not implemented in implicit_step.'
         stop
      end if

C     Collocation points
      dx = xout - xin
      do i=1,nu
         x(i) = xin + dx * c(i)
      end do

C     Stages
      do i=1,nu
         do j=1,ny
            xi(j,i) = yin(j)
            do ii=1,i-1
               xi(j,i) = xi(j,i) + dx * a(i,ii) * f(j,ii)
            end do
         end do
C        Evaluate derivatives f at stage xi
C        Inner patch: junk is irrelevant
C        Outer patch: junkc=Delta, junkf1=xi, junkf2=dxidtau
C        Future patch: junk
         call derivs(ny, xi(1,i), x(i), f(1,i), junkc, junkf1, junkf2)
      end do

C     Evaluate new value of y
      do j=1,ny
         yout(j) = yin(j)
         do i=1,nu
            yout(j) = yout(j) + dx * b(i) * f(j,i)
         end do
      end do
C     There is no iterative process here
      itsreach = 1

      return

      end

