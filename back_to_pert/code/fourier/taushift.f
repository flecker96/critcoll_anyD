
      program taushift

      implicit none
      integer nymax, ny, i
      include '../nymax.inc'
      double precision xin(nymax), xout(nymax), caux(2*nymax),
     $     a, b, phase, tau0, Delta, tau, tau0test

      ny = 64
      Delta = 1.d0

      open(unit=10,file='fc.dat',status='old')
      do i=1,ny
         read(10,*) xin(i)
      end do
      close(10)

      do i=1,ny
         caux(2*i-1) = xin(i) / ny
         caux(2*i) = 0.d0
      end do
      call four1(caux, ny, 1)
      a = caux(5)
      b = caux(6)
      phase = asin(a / sqrt(a**2 + b**2))
      tau0 = - phase * Delta / 2.d0 / 6.28318530717959d0
      write(6,99) tau0

      do i=1,ny
         tau = Delta * dble(i-1) / dble(ny) + tau0
         call fourinterpol(ny, xin, Delta, tau, xout(i))
      end do
      open(unit=10,file='fc.out',status='unknown')
      do i=1,ny
         write(10,99) xout(i)
      end do
      close(10)

      do i=1,ny
         caux(2*i-1) = xout(i) / ny
         caux(2*i) = 0.d0
      end do
      call four1(caux, ny, 1)
      a = caux(5)
      b = caux(6)
      phase = asin(a / sqrt(a**2 + b**2))
      tau0test = - phase * Delta / 2.d0 / 6.28318530717959d0
      write(6,99) tau0test

      open(unit=10,file='Up.dat',status='old')
      do i=1,ny
         read(10,*) xin(i)
      end do
      close(10)
      do i=1,ny
         tau = Delta * dble(i-1) / dble(ny) + tau0
         call fourinterpol(ny, xin, Delta, tau, xout(i))
      end do
      open(unit=10,file='Up.out',status='unknown')
      do i=1,ny
         write(10,99) xout(i)
      end do
      close(10)

      open(unit=10,file='psic.dat',status='old')
      do i=1,ny
         read(10,*) xin(i)
      end do
      close(10)
      do i=1,ny
         tau = Delta * dble(i-1) / dble(ny) + tau0
         call fourinterpol(ny, xin, Delta, tau, xout(i))
      end do
      open(unit=10,file='psic.out',status='unknown')
      do i=1,ny
         write(10,99) xout(i)
      end do
      close(10)

99    format(F24.16)

      return
      end
