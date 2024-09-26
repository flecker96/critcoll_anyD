
      program changeresolution

      implicit none
      integer nymax, nyin, nyout, i
      include '../nymax.inc'
      double precision x, caux(2*nymax)

      open(unit=10,file='changeresolution.par',status='old')
      read(10,*) nyin
      read(10,*) nyout
      close(10) 

      if (nyin .gt. nymax .or. nyout .gt. nymax) stop 'nymax too small'

      do i=1,nyin
         read(5,*) x
         caux(2*i-1) = x / nyin
         caux(2*i) = 0.d0
      end do

      call four1(caux, nyin, 1)

      if (nyin .eq. 2*nyout) then
         call halve(caux, caux, nyout)
      else if (nyout .eq. 2*nyin) then
         call double(caux, caux, nyin)
      else
         stop 'can only double or halve resolution'
      end if

      call four1(caux, nyout, -1)

      do i=1,nyout
         write(6,99) caux(2*i-1)
      end do

99    format(F24.16)

      return
      end
