C***********************************************************************
C
C     Program for writing the Choptuik spacetime
C     from converged initial data to files. 
C
C     INPUT:
C
C        "Delta.out", "fc.out", "psic.out", "Up.out"
C
C     OUTPUT:
C
C        The spacetime between xcut0 and xcut1 stored in
C
C        fB.dat, uB.dat, vB.dat, ia2B.dat
C
C        taylor_coeff_R.dat
C
C        delfc.dat  , delpsic.dat,  delUp.dat 
C
C        The latter will be used as initial data for the perturbations.
C
C***********************************************************************


      program bg_to_file

C     *******************************
C     **** Variable declarations ****
C     *******************************

      implicit none
      
C     tau related variables
      integer ny, n3, nymax, n3max, j, tstep
      include '../nymax.inc'
      parameter(n3max = 3*nymax/4)
      double precision in0(n3max), out0(n3max),
     $     Delta, fc(nymax), Up(nymax), psic(nymax),
     $     UpB(nymax), psicB(nymax), fcB(nymax)
 
C     x related variables
      integer nx, im, ileft, imid, iright, i, outevery,
     $     nxmax, itsreach, nleft, nright
      parameter(nxmax = 100001)
      double precision xc, xm, xp, xleft, xmid, xright,
     $     prec_irk, xxp(nxmax), dx,
     $     zleft, zm, zright, dz, zz, tol

C     Newton method variables
      integer maxits
      parameter(maxits=100)
      double precision eps_newt, d, err, prec_newt, slowerr, fac

C     Output
      logical verbose, debug, useloggrid, gridexists

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
      read(10,*) useloggrid
      read(10,*) nleft
      read(10,*) nright
      read(10,*) tol
      read(10,*) tstep
      read(10,*) prec_irk
      read(10,*) debug
      close(10)
     
C     Create directory where the files necessary for the lincode
C     are stored, remove if already present
      call system('rm -rf bg_data; mkdir bg_data')
     
C     ny is already the doubled number of modes (anti aliasing)
      n3 = 3*ny/4

C     Limits of the inner patch
      xc = 0.d0
      xp = 1.d0
C     Midpoint. 
      xm = xmid

C     **************
C     **** Grid ****
C     **************

C     Check for already existing grid.
      inquire(file='grid.dat', exist=gridexists)

      if (gridexists) then 
C        Load already existing grid
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

C        Copy grid files to bg_data diretory
         call system('cp grid.par bg_data/grid.par')
         call system('cp grid.dat bg_data/grid.dat')

      else if (useloggrid) then

C       Grid points:
C       ileft: leftmost point of the grid, corresponding to xleft
C       im: change of logarithmic axes, corresponding to xm
C       imid: junction point between left and right evolutions
C       iright: rightmost point of the grid, corresponding to xright

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
         open(unit=10,file='bg_data/loggrid.dat',status='new')
         do i=1,nx
            if (i.ne.nx) then
               dx = xxp(i+1) - xxp(i)
            else
               dx = xxp(nx) - xxp(nx-1)
            end if
            write(10,'(G23.16,X,G23.16)') xxp(i), dx
         end do
         close(10)
         write(6,*) 'INFO: Grid: ileft, imid, iright: ',
     $    ileft, imid, iright
         write(6,*) 'INFO: Grid: xleft, xmid, xright: ', 
     $    xxp(ileft), xxp(imid), xxp(iright)
     
      else
         stop 'No grid.dat files found.'
      end if

C     *******************
C     ****** Data *******
C     *******************

      write(6,*) 'INFO: Starting from out files'
C     Read in free functions from out files of converged BG code. 
      open(unit=10, file='Delta.out', status='old')
      read(10,*) Delta
      close(10)
      open(unit=10, file='fc.out', status='old')
      do j=1,ny
         read(10,*) fc(j)
      end do
      close(10)
      open(unit=10, file='psic.out', status='old')
      do j=1,ny
         read(10,*) psic(j)
      end do
      close(10)    
      open(unit=10, file='Up.out', status='old')
      do j=1,ny
         read(10,*) Up(j)
      end do
      close(10)
       
C     Write Delta to file         
      open(unit=10, file='bg_data/Delta.dat', status='new')
            write(10,99) Delta
      close(10)
         
C     Distinguish between downsampling and no downsampling in tau
      if (tstep.eq.1) then
         open(unit=10, file='bg_data/fcB.dat', status='new')
         do j=1,ny
            write(10,99) fc(j)
         end do
         close(10)
         open(unit=10, file='bg_data/psicB.dat', status='new')
         do j=1,ny
            write(10,99) psic(j)
         end do
         close(10)
         open(unit=10, file='bg_data/UpB.dat', status='new')
         do j=1,ny
            write(10,99) Up(j)
         end do
         close(10)
      else if (tstep.eq.2) then
         open(unit=10, file='bg_data/fcB.dat', status='new')
         call changeresolution(ny,ny/2,fc,fcB)
         do j=1,ny/2
            write(10,99) fcB(j)
         end do
         close(10)
         open(unit=10, file='bg_data/psicB.dat', status='new')
         call changeresolution(ny,ny/2,psic,psicB)
         do j=1,ny/2
            write(10,99) psicB(j)
         end do
         close(10)
         open(unit=10, file='bg_data/UpB.dat', status='new')
         call changeresolution(ny,ny/2,Up,UpB)
         do j=1,ny/2
            write(10,99) UpB(j)
         end do
         close(10)
      else 
         write(6,*) 'Can only downsample by factor of 2.'
         stop
      end if
                  
C     Use a scaled version of converged BG data as IVs for perturbations
      if (tstep.eq.1) then
         open(unit=10, file='bg_data/delfc.dat', status='new')
         do j=1,ny
            write(10,99) fc(j) * 1.d-3
         end do
         close(10)
         open(unit=10, file='bg_data/delpsic.dat', status='new')
         do j=1,ny
            write(10,99) psic(j) * 1.d-3
         end do
         close(10)
         open(unit=10, file='bg_data/delUp.dat', status='new')
         do j=1,ny
            write(10,99) Up(j) * 1.d-3
         end do
         close(10)
      else if (tstep.eq.2) then
         open(unit=10, file='bg_data/delfc.dat', status='new')
         do j=1,ny/2
            write(10,99) fcB(j) * 1.d-3
         end do
         close(10)
         open(unit=10, file='bg_data/delpsic.dat', status='new')
         do j=1,ny/2
            write(10,99) psicB(j) * 1.d-3
         end do
         close(10)
         open(unit=10, file='bg_data/delUp.dat', status='new')
         do j=1,ny/2
            write(10,99) UpB(j) * 1.d-3
         end do
         close(10)
      else 
         write(6,*) 'Can only downsample by factor of 2.'
         stop
      end if

      
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

    
      write(6,*) 'INFO: Nt of output: ',ny/(2*tstep)          
      

C     Force log output
      call flush(6)

C     ******************
C     **** Shooting ****
C     ******************
            
      write(6,*) 'INFO: Shooting once and extracting fields... '
      call shoot_to_file(ny, nx, d, ileft, iright, imid,
     $      n3, in0, out0, xxp, tstep, prec_irk, debug, 
     $      itsreach)
      
C     Calculate new error
      err = 0.d0
      do j=1,n3
         err = err + out0(j)**2
      end do
      err = sqrt(err/n3)
      fac = min(1.d0,slowerr/err)
      write(*,*) 'Mismatch: ', err
      write(*,*) 'Done.'

99    format(F24.16)
         
      end program

      
