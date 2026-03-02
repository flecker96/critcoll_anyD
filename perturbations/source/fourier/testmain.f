      program testmain

      implicit none
      integer ny, j
      parameter (ny=64)
      double precision Delta, Up(ny), dUpdtau(ny)

C     Read in initial data at the past light cone. 
      open(unit=10, file='Delta.dat', status='old')
      read(10,*) Delta
      close(10)
      open(unit=10, file='Up.dat', status='old')
      do j=1,ny
         read(10,*) Up(j)
      end do
      close(10)

      call fourdiff1(ny, Up, dUpdtau, Delta)

      open(unit=10, file='dUpdtau.junk', status='unknown')
      do j=1,ny
         write(10,*) dUpdtau(j)
      end do
      close(10)

      end
