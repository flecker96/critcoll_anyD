
      subroutine changeresolution(nyin, nyout, fin, fout)

      implicit none
      integer nymax, nyin, nyout, i
      include '../nymax.inc'
      double precision fin(nyin), fout(nyout), caux(2*nymax)

      
      if (nyin .gt. nymax .or. nyout .gt. nymax) stop 'nymax too small'

      do i=1,nyin
         caux(2*i-1) = fin(i) / nyin
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
         fout(i) = caux(2*i-1)
      end do

      return
      end
