C***********************************************************************
C
C     Main program for calculating the Choptuik spacetime
C     (inner patch), modified version of code by CG & JM 
C
C     Structure of the program
C
C        1.- Variable declarations
C        2.- Read parameters
C        3.- Input free data or previous checkpointing
C        4.- Generate grid
C        5.- Shooting and Newton's method
C
C     Every (non tmp) output is controlled by the output subroutine.
C
C     Units:
C        tmp files: 10
C        output: 11-20
C        debug leftdata: 30-38, 40-48, 51-54
C        debug rightdata: 60-73
C        shoot1: 79-83
C        more debug: 90-99
C
C***********************************************************************

      program shootmain_inner

C     *******************************
C     **** Variable declarations ****
C     *******************************

      implicit none
      
C     tau related variables
      integer ny, n3, nymax, n3max, j
      include '../nymax.inc'
      parameter(n3max = 3*nymax/4)
      double precision in0(n3max), out0(n3max), change(n3max),
     $     Delta, fc(nymax), Up(nymax), psic(nymax), 
     $     outjunk(n3max)
 
C     x related variables
      integer nx, im, ileft, imid, iright, i, outevery, method,
     $     nxmax, itsreach, nleft, nright
      parameter(nxmax = 30001)
      double precision xc, xm, xp, xleft, xmid, xright,
     $     prec_irk, xxp(nxmax),
     $     zleft, zm, zright, dz, zz, tol, dxini

C     Newton method variables
      integer its, indx(n3max), ivar, maxits
      parameter(maxits=100)
      double precision in1(n3max), out1(n3max),
     $     doutdin(n3max,n3max), doutdinLU(n3max,n3max),
     $     eps_newt, d, err, prec_newt, slowerr, fac, b, errold

C     Output
      integer itsread, ivarread
      logical verbose, debug, CPexists, gridexists, useloggrid

C     ********************
C     **** Parameters ****
C     ********************

C     Read parameters from file
      open(unit=10,file='shoot_inner.par',status='old')
      read(10,*) ny
      read(10,*) d
      read(10,*) xleft
      read(10,*) xmid
      read(10,*) xright
      read(10,*) eps_newt
      read(10,*) prec_newt
      read(10,*) verbose
      read(10,*) slowerr
      read(10,*) outevery
      read(10,*) nleft    !Only needed for GG grid
      read(10,*) nright   !Only needed for GG grid
      read(10,*) tol  
      read(10,*) method  !Not needed --> tstep
      read(10,*) prec_irk
      read(10,*) debug
      close(10)

C     Parameters for adaptive stepsize algorithm
      dxini = 1.d-6
      useloggrid = .false.

C     ny is already the doubled number of modes (anti aliasing)
C     n3 is the number of variables in the outer Newton alg
      n3 = 3*ny/4

C     Limits of the inner patch
      xc = 0.d0
      xp = 1.d0
C     Midpoint. Change of logarithmic axes.
      xm = xmid
     
C     *******************
C     **** Free data ****
C     *******************

C     Check for checkpointing of previous variations.
      inquire(file='status.junk', exist=CPexists)
      if (CPexists) then

C     Prepare free data from previous run.
         write(6,*) 'INFO: Previous variations found'

C        Scan status.junk to find out how far we have already got.
         open(unit=10, file='status.junk', status='old')
            do its=1,maxits
               do ivar=1,n3
                  read(10,*, end=88) itsread, ivarread
                  if ((its.eq.itsread).and.(ivar.eq.ivarread)) then
                     continue
                  else
                     write(6,*) 'ERROR: Inconsistent status.junk file:'
                  end if
               end do
            end do
            write(6,*) 'ERROR: EOF not found'
C        Code only branches here when EOF of status.junk is found.
 88      close(10)

C        Decide whether we have full iterations or not.
         if (ivarread.eq.n3) then
C           Yes. Full iterations completed. Do nothing.
            write(6,*) 'INFO: Found ', itsread, ' full iterations.'
         else
C           No. Incomplete iteration found. Read it.
            itsread = itsread - 1
            write(6,*) 'INFO: Found ', itsread, ' iterations and ',
     $         ivarread, ' variations.'
            if (itsread.eq.-1) then
               write(6,*) 'ERROR: First iteration is incomplete.'
               write(6,*) 'ERROR: Please remove and restart.'
               stop
            end if
C           Read in doutdin.
            open(unit=10, file='doutdin.junk', status='old')
C           Skip previous full iterations.
            do its=1,itsread
               do ivar=1,n3
                  do j=1,n3
                     read(10,*)
                  end do
               end do
            end do
C           Read already computed variations from the current iteration.
            do ivar=1,ivarread
               do j=1,n3
                  read(10,*) doutdin(j,ivar)
                  doutdinLU(j,ivar) = doutdin(j,ivar)
               end do
            end do 
            close(10)
         end if

C        Read in in0.
         open(unit=10, file='in.junk', status='old')
         do its=1,itsread+1
            do j=1,n3
               read(10,*) in0(j)
            end do
         end do
         close(10)

      else

C     Brand new evolution.
         itsread = 0
         ivarread = 0
         write(6,*) 'INFO: Starting from dat files'
C     Read in free functions from dat files. 
         open(unit=10, file='Delta.dat', status='old')
         read(10,*) Delta
         close(10)
         open(unit=10, file='fc.dat', status='old')
         do j=1,ny
            read(10,*) fc(j)
         end do
         close(10)
         open(unit=10, file='psic.dat', status='old')
         do j=1,ny
            read(10,*) psic(j)
         end do
         close(10)    
         open(unit=10, file='Up.dat', status='old')
         do j=1,ny
            read(10,*) Up(j)
         end do
         close(10)
         
C     Pack them into in0: [ FT(psic), FT(Up), FT(fc) ], where each of
C     the Fourier transforms has n3/3=ny/4 elements and contains the
C     upper half of the Fourier spectrum of the variables. (Thus, ny
C     must be a multiple of 8.) Odd functions are trivially defined.
C     Even functions have 0 imaginary part in the constant mode, which
C     is replaced by the amplitude of the high-frequency cosine.
         call mypack(n3, in0, ny, Up, psic, fc)
         
C     We fix a global phase assuming that fc has no fundamental cosine 
C     at x=xc. (Re(f2)=0) This mode is translated into in0(ny/2+3). Replace it by
C     Delta to make clear that we have n3 free variables.
         write(6,*) 'INFO: Low frequency cosine of fc:', in0(ny/2+3)
         in0(ny/2+3) = Delta

      end if

      write(6,'(A11,I4,A6,F6.4)') ' INFO: Nt = ',ny/2,', d = ',d

C     Force log output
      call flush(6)


C     ******************
C     ***** Grid *******
C     ******************

C     Check for already existing grid.
      inquire(file='grid.dat', exist=gridexists)

      if (useloggrid) then
C        This is only used if we want to have the old definition
C        of the grid according to the GG paper

C        Point index
         ileft = 1
         im = nleft + 1
         iright = nleft + nright + 1

C        Log grid
         zleft  = log(xleft -xc) - log(xp-xleft)
         zm = log(xmid-xc) - log(xp-xmid)
         zright = log(xright-xc) - log(xp-xright)
      
C        Set grid
         nx = iright
         if (nx.gt.nxmax) stop 'nxmax too small in shootmain'

         xxp(ileft) = xleft

         dz = (zm - zleft) / dble(nleft)
         do i=ileft+1,im-1
            zz = zleft + dble(i-ileft) * dz
            xxp(i) = (xc + xp * exp(zz) )/ ( 1.d0 + exp(zz) )
         end do

         xxp(im) = xm

         dz = (zright - zm) / dble(nright)
         do i=im+1,iright-1
            zz = zm + dble(i-im) * dz
            xxp(i) = (xc + xp * exp(zz) )/ ( 1.d0 + exp(zz) )
         end do

         xxp(iright) = xright
      
C        Look for imid
         if ((xmid.lt.xleft).or.(xmid.gt.xright)) then
            write(6,*) 'ERROR: xmid is out of range.'
            stop
         end if
         if (xmid.eq.xleft) imid = 1
         do i=ileft+1,iright
            if ((xxp(i).ge.xmid).and.(xxp(i-1).lt.xmid)) then
               xmid = xxp(i)
               imid = i
            end if
         end do

C        Output grid
         open(unit=10,file='loggrid.dat',status='unknown')
         do i=1,nx
            write(10,'(G23.16)') xxp(i)
         end do
         close(10)
         write(6,*) 'INFO: Grid: ileft, imid, iright: ',
     $      ileft, imid, iright
         write(6,*) 'INFO: Grid: xleft, xmid, xright: ', 
     $      xxp(ileft), xxp(imid), xxp(iright)

      else if (gridexists) then
C     Load already existing grid
         write(6,*) 'INFO: Found grid.dat file, loading grid...'

         open(unit=10,file='grid.dat',status='old')
         do i=1,nxmax
            read(10,'(G23.16)',end=89) xxp(i)
         end do
 89      close(10)

         open(unit=10,file='grid.par',status='old')
            read(10,'(I5,G23.16)') ileft, xxp(ileft)
            read(10,'(I5,G23.16)') imid, xxp(imid)
            read(10,'(I5,G23.16)') iright, xxp(iright)
         close(10)

C        Point index
         nx = iright

         write(6,*) 'INFO: Grid: ileft, imid, iright: ',
     $      ileft, imid, iright
         write(6,*) 'INFO: Grid: xleft, xmid, xright: ', 
     $      xxp(ileft), xxp(imid), xxp(iright)
      else
C        Generate a new grid with adapive stepsize algorithm
C        Uses xleft, xright, xmid and input and gives
C        xxp, nleft, nright as output
         write(6,90) ' INFO: Generating new x-grid with LTE = ', tol
 90      format (A,ES15.2)    

         call gen_grid(ny, d, xleft, xright, xmid,
     $        n3, in0, outjunk, prec_irk, tol, xxp, nleft, nright, 
     $        dxini, debug)

C        Point index
         ileft = 1
         imid = nleft + 1
         iright = nleft + nright + 1
         nx = iright
         
C        Write grid to file
         open(unit=10,file='grid.dat',status='new')
         do i=1,nx
            write(10,'(G23.16)') xxp(i)
         end do
         close(10)
         open(unit=10,file='grid.par',status='new')
            write(10,'(I5,G23.16)') ileft, xxp(ileft)
            write(10,'(I5,G23.16)') imid, xxp(imid)
            write(10,'(I5,G23.16)') iright, xxp(iright)
         close(10)
         write(6,*) 'INFO: Grid: ileft, imid, iright: ',
     $      ileft, imid, iright
         write(6,*) 'INFO: Grid: xleft, xmid, xright: ', 
     $      xxp(ileft), xxp(imid), xxp(iright)

      end if


C     ******************
C     **** Shooting ****
C     ******************

      if (verbose)
     $  call output('open', ny, d, n3, in0, out0, change, err, fac)

C     Shooting
      itsreach = 1

C     Newton iterations
      do its=itsread,maxits

         if (its.ge.(itsread+1)) then
            errold = err
         else
            errold = 1.d0
         end if

         write(6,*) 'INFO: Shooting iteration ', its
         call shoot1(ny, nx, d, ileft, iright, imid, 0,
     $      n3, in0, out0, outevery, xxp, prec_irk, 
     $      debug, itsreach)
         
C        Calculate new error
         err = 0.d0
         do j=1,n3
            err = err + out0(j)**2
         end do
         err = sqrt(err/n3)
         fac = min(1.d0,slowerr/err)
         write(*,*) 'Mismatch: ', err
         if (verbose) then
            if ((its.eq.itsread).and.CPexists) then
C              Do not ouptut, to avoid duplicating initial output.
            else
               call output('write', ny, d, n3, in0 
     $                      ,out0, change, err, fac)
            end if
         end if

C        Debug stop
         if (debug) stop
         
C        If error is small enough, branch out
         if (err .lt. prec_newt) goto 1
        
C        If mismatch does not change any more, branch out as well 
C         if (INT(LOG10(err)).eq.INT(LOG10(errold))) goto 1     
         
C        Else, perform new iteration.
         do ivar = 1,n3
C         do ivar=ivarread+1,n3

C           Store input free data
            do j=1,n3
               in1(j) = in0(j)
            end do
C           Change variable ivar
            in1(ivar) = in1(ivar) + eps_newt
            
C           Output info
            write(6,*) '  Max its reached:', itsreach, '. Varying', ivar
            call flush(6)

C           Shoot new free data
            itsreach = 1
            call shoot1(ny, nx, d, ileft, iright, imid, ivar,
     $           n3, in1, out1, 0, xxp, prec_irk, 
     $           debug, itsreach)

C           Calculate derivative with respect to variable ivar
            do j=1,n3
               doutdin(j,ivar) = (out1(j) - out0(j)) / eps_newt
               doutdinLU(j,ivar) = doutdin(j,ivar)
            end do
            if (verbose) 
     $         call output_step(n3, its+1, ivar, doutdin(1,ivar))

         end do
         
C        Calculate suggested change in doutdin.change = out0
C        LU decomposition of doutdin, stored in the same space
         call ludcmp(doutdinLU, n3, n3max, indx, b)
C        Feed out0 to the subroutine using change
         do j=1,n3
            change(j) = out0(j)
         end do
         call lubksb(doutdinLU, n3, n3max, indx, change)
C        Improve result
         call mprove(doutdin, doutdinLU, n3, n3max, indx, out0, change)
         call mprove(doutdin, doutdinLU, n3, n3max, indx, out0, change)
         call mprove(doutdin, doutdinLU, n3, n3max, indx, out0, change)
C        New free data
         do j=1,n3
            in0(j) = in0(j) - fac * change(j)
         end do

C        Iterate Newton method
         ivarread = 0

      end do

C     The code only branches here once err < prec_newt.
 1    continue

C     Final output
      if (verbose) then
         call output('write_fin', ny, d, n3, in0, out0
     $                       , change, err, fac)
         call output('close', ny, d, n3, in0, out0
     $                , change, err, fac)
      end if

      end



C***********************************************************************


      subroutine output(action, ny, d, n3, in0, out0, change, err, fac)

      implicit none

      character action*(*)
      integer ny, n3
      double precision in0(n3), out0(n3), change(n3), err, fac, d

      integer nymax, j
      include '../nymax.inc'
      double precision Delta, Up(nymax), fc(nymax), psic(nymax),
     $   y(nymax), vp(nymax), junk(nymax), xp
      logical CPexists

      if (action.eq.'open') then
         inquire(file='status.junk', exist=CPexists)
         if (CPexists) then 
            open(unit=11, file='Delta.junk', 
     $            position='append', status='old')
            open(unit=12, file='in.junk', 
     $            position='append', status='old')
            open(unit=13, file='out.junk',
     $            position='append', status='old')
            open(unit=14, position='append', 
     $            file='doutdin.junk', status='old')
            open(unit=15, file='change.junk', 
     $            position='append', status='old')
            open(unit=16, file='err.junk', 
     $            position='append', status='old')
            open(unit=17, file='fc.junk', 
     $            position='append', status='old')
            open(unit=18, file='psic.junk', 
     $            position='append', status='old')
            open(unit=19, file='Up.junk', 
     $            position='append', status='old')
            open(unit=20, position='append', 
     $            file='status.junk', status='old')


         else
            open(unit=11, file='Delta.junk', status='new')
            open(unit=12, file='in.junk', status='new')
            open(unit=13, file='out.junk', status='new')
            open(unit=14, file='doutdin.junk', status='new')
            open(unit=15, file='change.junk', status='new')
            open(unit=16, file='err.junk', status='new')
            open(unit=17, file='fc.junk', status='new')
            open(unit=18, file='psic.junk', status='new')
            open(unit=19, file='Up.junk', status='new')
            open(unit=20, file='status.junk', status='new')
         end if
      else if (action.eq.'write') then

         Delta = in0(ny/2+3)
         in0(ny/2+3) = 0.d0
         call myunpack(n3, in0, ny, Up, psic, fc)
         in0(ny/2+3) = Delta

         write(11,99) Delta
         call flush(11)

         do j=1,n3
            write(12,99) in0(j)
         end do
         write(12,*) 
         call flush(12)

         do j=1,n3
            write(13,99) out0(j)
         end do
         write(13,*) 
         call flush(13)

         do j=1,n3
            write(15,99) change(j)
         end do
         write(15,*) 
         call flush(15)

         write(16,*) err, fac
         call flush(16)

         do j=1,ny
            write(17,99) fc(j)
         end do
         write(17,*)
         call flush(17)

         do j=1,ny
            write(18,99) psic(j)
         end do
         write(18,*)
         call flush(18)

         do j=1,ny
            write(19,99) Up(j)
         end do
         write(19,*)
         call flush(19)

      else if (action.eq.'write_fin') then
         Delta = in0(ny/2+3)
         in0(ny/2+3) = 0.d0
         call myunpack(n3, in0, ny, Up, psic, fc) 
         in0(ny/2+3) = Delta

         open(unit=10, file='Delta.out', status='unknown')
         write(10,99) Delta
         close(10)

         open(unit=10, file='fc.out', status='unknown')
         do j=1,ny
            write(10,99) fc(j)
         end do
         close(10)

         open(unit=10, file='psic.out', status='unknown')
         do j=1,ny
            write(10,99) psic(j)
         end do
         close(10)     

         open(unit=10, file='Up.out', status='unknown')
         do j=1,ny
            write(10,99) Up(j)
         end do
         close(10)

C     Calculate Vp from Delta and Up.
         xp = 1.d0
         call rightdata(ny, d, xp, Delta, Up, y, .FALSE.)
         call fieldsfromy(ny, d, y, xp,
     $        junk, vp, junk, junk, Delta,
     $        junk, junk, junk) 

         open(unit=10, file='Vp.out', status='unknown')
         do j=1,ny
            write(10,99) vp(j)
         end do
         close(10)

      else if (action.eq.'close') then
         do j=11,20
            close(j)
         end do
         write(*,*) 'INFO: Files closed. Evolution ended.'
      else
         write(*,*) 'ERROR: Unrecognized action ', action, 'in output.'
      end if

 99   format(G24.16)

      return
      end

C***********************************************************************


      subroutine output_step(n3, its, ivar, doutdin)

      implicit none

      integer n3, its, ivar
      double precision doutdin(n3)

      integer j

C     Derivative of variable j of out0 with respect to variable ivar
C     of in0
  
      do j=1,n3
         write(14,99) doutdin(j)
      end do
      write(14,*) 
      call flush(14)
      write(20,*) its, ivar
      call flush(20)

 99   format(G24.16)

      return
      end

      
