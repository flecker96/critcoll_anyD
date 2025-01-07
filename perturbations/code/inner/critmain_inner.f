C***********************************************************************
C
C     Main program for calculating the critical exponent for the
C     Choptuik spacetime (inner patch) 
C
C     Structure of the program
C
C        1.- Variable declarations
C        2.- Read parameters
C        3.- Read grid and background data
C        4.- Shooting and root finding
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


      program crit_inner

C     *******************************
C     **** Variable declarations ****
C     *******************************

      implicit none
      
C     tau related variables
      integer ny, n3, nymax, n3max, j, tstep
      include '../nymax.inc'
      parameter(n3max = 3*nymax/4)
      double precision in0(n3max),
     $      fc(nymax), Up(nymax), psic(nymax), gam
 
C     x related variables
      integer nx, ileft, imid, iright, i, outevery,
     $     nxmax, nleft, nright
      parameter(nxmax = 20001)
      double precision xc, xp, xleft, xmid, xright,
     $     prec_irk, xxp(nxmax), dx

C     Background variables
      double precision Delta,  
     $     uB(nymax, nxmax), vB(nymax,nxmax), 
     $     fB(nymax,nxmax), ia2B(nymax, nxmax),
     $     fcB(nymax), psicB(nymax), UpB(nymax), in0B(n3max)

C     Newton method variables
      integer its, ivar, maxits
      parameter(maxits=100)
      double precision 
     $     doutdin(n3max,n3max),
     $     eps, d, errmax, slowerr,
     $     gamlist(maxits), detlist(maxits),
     $     gamstep, det, fac

C     Brent method variables
      double precision x1, x2, x3, dx1, dx2,
     $   f1, f2, f3, tolbrent, p, q, r, s, xm,
     $   brenttol

C     Output
      integer itsread, ivarread
      logical verbose, debug, CPexists, foundzero, newfac, gridexists

C     ********************
C     **** Parameters ****
C     ********************

C     Read parameters from file
      open(unit=10,file='perturbations.par',status='old')
      read(10,*) ny
      read(10,*) d
      read(10,*) xleft
      read(10,*) xmid
      read(10,*) xright
      read(10,*) eps
      read(10,*) errmax
      read(10,*) verbose
      read(10,*) slowerr
      read(10,*) outevery
      read(10,*) nleft
      read(10,*) nright
      read(10,*) tstep
      read(10,*) prec_irk
      read(10,*) debug
      read(10,*) brenttol
      read(10,*) gamstep
      close(10)
     
C     ny is already the doubled number of modes (anti aliasing)
C     Here it is additionally divided by tstep (takes values 1 or 2)
      ny = ny / tstep
      n3 = 3*ny/4
      
C     Limits of the inner patch
      xc = 0.d0
      xp = 1.d0

C     **************
C     **** Grid ****
C     **************


C     Check for already existing grid.
      inquire(file='bg_data/grid.dat', exist=gridexists)

      if (gridexists) then 
C        Load already existing grid
         write(6,*) 'INFO: Found grid.dat file, loading grid...'

         open(unit=10,file='bg_data/grid.dat',status='old')
         do i=1,nxmax
            read(10,'(G23.16)',end=89) xxp(i)
         end do
 89      close(10)

         open(unit=10,file='bg_data/grid.par',status='old')
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
C        Check if loggrid is there         
         inquire(file='bg_data/loggrid.dat', exist=gridexists)
         if (.not.gridexists) stop 'No grid file found.'
         write(6,*) 'INFO: Loading loggrid.dat file...'

C        Point index
         ileft = 1
         iright = nleft + nright + 1

         nx = iright

C        Read in x-grid
         open(unit=10,file='bg_data/loggrid.dat',status='old')
            do i=1,nx
               read(10,*) xxp(i), dx
            end do
         close(10)
      
C        Look for imid
         do i=ileft+1,iright
            if ((xxp(i).ge.xmid).and.(xxp(i-1).lt.xmid)) then
               xmid = xxp(i)
               imid = i
            end if
         end do

         write(6,*) 'INFO: Grid: ileft, imid, iright: ',
     $    ileft, imid, iright
         write(6,*) 'INFO: Grid: xleft, xmid, xright: ', 
     $    xxp(ileft), xxp(imid), xxp(iright)
     
      end if

C     *******************   
C     *** Background ****
C     *******************

C     Read in fields
      write(6,*) 'INFO: Reading in background'
      open(unit=10,file='bg_data/uB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) uB(j,i)
            end do
         end do
      close(10)
      
      open(unit=10,file='bg_data/vB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) vB(j,i)
            end do
         end do
      close(10)
      
      open(unit=10,file='bg_data/fB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) fB(j,i)
            end do
         end do
      close(10)
      
      open(unit=10,file='bg_data/ia2B.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) ia2B(j,i)
            end do
         end do
      close(10)
         

C     Read in background initial data from out files. 
      open(unit=10, file='bg_data/Delta.dat', status='old')
        read(10,*) Delta
      close(10)
      open(unit=10, file='bg_data/fcB.dat', status='old')
        do j=1,ny
          read(10,*) fcB(j)
        end do
      close(10)
      open(unit=10, file='bg_data/psicB.dat', status='old')
        do j=1,ny
          read(10,*) psicB(j)
        end do
      close(10)    
      open(unit=10, file='bg_data/UpB.dat', status='old')
        do j=1,ny
           read(10,*) UpB(j)
        end do
      close(10)
       
      
C     ***************************
C     **** Perturbation data ****
C     ***************************

C     Read in free functions from dat files. 
      open(unit=10, file='gamma.dat', status='old')
      read(10,*) gam
      close(10)
      open(unit=10, file='bg_data/delfc.dat', status='old')
      do j=1,ny
         read(10,*) fc(j)
      end do
      close(10)
      open(unit=10, file='bg_data/delpsic.dat', status='old')
      do j=1,ny
         read(10,*) psic(j)
      end do
      close(10)    
      open(unit=10, file='bg_data/delUp.dat', status='old')
      do j=1,ny
         read(10,*) Up(j)
      end do
      close(10)
      
C     Initialization for root finding methods
      do its=1,maxits
         gamlist(its) = 0.d0
         detlist(its) = 0.d0
      end do
      
      gamlist(1) = gam     

      fac = 1.d0
      foundzero = .false.
      
C     This removes the old iderations; makes the part with 
C     CPexists = .true. obsolete because that part is
C     not yet compatible with the Brent algorithm 
C     -> possible future implementation
      call system('rm -rf pert_junk')
      
C     Check for checkpointing of previous variations.
      inquire(file='pert_junk/status.junk', exist=CPexists)
      if (CPexists) then

C     Prepare free data from previous run.
         write(6,*) 'INFO: Previous variations found'

C        Scan status.junk to find out how far we have already got.
         open(unit=10, file='pert_junk/status.junk', status='old')
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
C        Read in doutdin.
            open(unit=10, file='pert_junk/doutdin.junk', status='old')
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
               end do
            end do 
            close(10)
         end if         


C        Read in lamb and det
         open(unit=10, file='pert_junk/detofgam.junk', status='old')
         do its=1,itsread
               read(10,*) gamlist(its), detlist(its), fac, foundzero
         end do
         close(10)

      else

C        Brand new evolution.
         itsread = 0
         ivarread = 0
         
      end if
      

   
C     Pack them into in0: [ FT(psic), FT(Up), FT(fc) ], where each of
C     the Fourier transforms has n3/3=ny/4 elements and contains the
C     upper half of the Fourier spectrum of the variables. (Thus, ny
C     must be a multiple of 8.) Odd functions are trivially defined.
C     Even functions have 0 imaginary part in the constant mode, which
C     is replaced by the amplitude of the high-frequency cosine.
      call mypack(n3, in0, ny, Up, psic, fc)
      call mypack(n3, in0B, ny, UpB, psicB, fcB)
         
C     We fix a global phase assuming that fcB has no fundamental cosine 
C     at x=xc. (Re(fB2)=0) This mode is translated into in0B(ny/2+3). Replace it by
C     Delta.
      write(6,*) 'INFO: Low frequency cosine of fcB:', in0B(ny/2+3)
      in0B(ny/2+3) = Delta

      write(6,'(A11,I4,A6,F6.4)') ' INFO: Nt = ',ny/2,', d = ',d

C     Force log output
      call flush(6)




C     **************************************
C     **** Finding root of det(A)(lamb) ****
C     **************************************

      if (verbose)
     $  call output('open', gam, det, fac, foundzero) 

C     First, perform fixed steps, until sign changes:
      if (.not.foundzero) then

         write(6,*) 'INFO: Stepping towards root ' 
         write(6,'(A13,F6.4,A11,F6.4)') '  gamstart = ', gam, 
     $          ', gamstep = ', gamstep

C        Every iteration, compute new normalization of determinant      
         newfac = .true.

      do its=itsread,maxits-1
         
         if (its.gt.0) then                        
               gamlist(its+1) = gamlist(its) + gamstep
         end if
         
         write(6,'(A17,G24.16)') 'INFO: gamma now = ', gamlist(its+1)

         call detofgam(ny, nx, d, gamlist(its+1), ileft, iright, imid,
     $      n3, in0, in0B, outevery, xxp, prec_irk, debug, 
     $      uB, vB, fB, ia2B, ivarread, its, det, 
     $      doutdin, eps, fac, newfac)

         write(6,'(A16,G24.16,A2,G24.16,A2)') '(gamma, det) = (', 
     $          gamlist(its+1), ' ,', det,' )'

         detlist(its+1) = det
         
C        Check if sign has changed
         if (its.gt.0) then
            if ((detlist(its+1) / detlist(its)) .lt. 0.0d0) then
            
               foundzero = .true.
               write(6,*) 'INFO: Found root, switching to Brent method.'
               
            end if
         end if

C        Write current step result to file
         if (verbose) then
               call output('write', gamlist(its+1), det, fac, foundzero)
         end if
         
         ivarread = 0
         
         if (foundzero) then
            itsread = its+1
            goto 1
         end if

      end do
      stop 'No root found in maxits iterations.'
      end if


C     Brent method for finding zero of det(A)(lamb)
C     Source: Numerical recipes book.
C     Needs two starting values on either side of root
 1    if (foundzero) then

      write(6,*) 'INFO: Launch Brent algorithm.'
C     Do not compute new normalization any more, important for Brent alg. to work
      newfac = .false.

C     Initialize Brent algorithm
      x1=gamlist(itsread)
      x2=gamlist(itsread-1)
      f1=detlist(itsread)
      f2=detlist(itsread-1)

      if((f1.gt.0..and.f2.gt.0.).or.(f1.lt.0..and.f2.lt.0.)) then
        stop 'Root not bracketed in Brent.'
      end if

      x3=x2
      f3=f2
      dx2=0.d0

C     Iterations
      do its=itsread,maxits-1

C        If x3 is on same side as x2, set x3=x1       
         if((f2.gt.0..and.f3.gt.0.).or.(f2.lt.0..and.f3.lt.0.)) then
            x3=x1
            f3=f1
            dx1=x2-x1
            dx2=dx1
         end if

C        Make sure that x2 is the better guess, switch with x1 if not the case
         if(abs(f3).lt.abs(f2)) then
            x1=x2
            x2=x3
            x3=x1
            f1=f2
            f2=f3
            f3=f1
         end if

C        xm is (signed) distance between x2 and midpoint
         tolbrent=6.d-15*abs(x2) + brenttol / 2.d0
         xm=(x3-x2) / 2.d0

C        Write current step result to file
         if (its.eq.itsread) then
C           Do not ouptut, to avoid duplicating initial output.
         else
            gamlist(its) = x2
            detlist(its) = f2

            write(6,'(A16,G24.16,A2,G24.16,A2)') '(gamma, det) = (', 
     $             x2, ' ,', f2,' )'

            call output('write', gamlist(its), detlist(its),
     $                    fac, foundzero)
         end if

C        Break here if tolerance is satisfied
         if(abs(xm).le.tolbrent .or. f2.eq.0.) goto 2

         if(abs(dx2).ge.tolbrent .and. abs(f1).gt.abs(f2)) then
C           Attempt inverse quadratic interpolation            
            s=f2/f1
            if(x1.eq.x3) then
                p=2.d0*xm*s
                q=1.d0 - s
            else
                q=f1/f3
                r=f2/f3
                p=s*(2.d0*xm*q*(q-r)-(x2-x1)*(r-1.d0))
                q=(q-1.d0)*(r-1.d0)*(s-1.d0)
            endif

C           Check whether in bounds.
            if(p.gt.0.) q = -q

            p=abs(p)
            
            if(2.d0*p .lt. 
     $          min(3.d0*xm*q-abs(tolbrent*q),abs(dx2*q))) then
C               Accept interpolation.
                write(6,*) 'This step: Use inverse quadr interpolation.'
                dx2=dx1
                dx1=p/q
            else
C               Interpolation failed, use bisection.
                write(6,*) 'This step: Interpolation failed, use', 
     $                       'bisection.'
                dx1=xm
                dx2=dx1
            endif
         else
C        Bounds decreasing too slowly, use bisection.
            write(6,*) 'This step: Bounds decreasing too slowly,',
     $       'use bisection.'
            dx1=xm
            dx2=dx1
         endif

C        Move last best guess to x1.
         x1=x2
         f1=f2

C        Evaluate new trial root.
         if(abs(dx1) .gt. tolbrent) then
            x2=x2+dx1
         else
            x2=x2+sign(tolbrent,xm)
         endif

C        Evaluate f(x2)
         call detofgam(ny, nx, d, x2, ileft, iright, imid,
     $      n3, in0, in0B, outevery, xxp, prec_irk, debug, 
     $      uB, vB, fB, ia2B, ivarread, its, f2, 
     $      doutdin, eps, fac, newfac)

         ivarread = 0

      end do

      stop 'Exceeded maxits in Brent.' 
         
      end if

C     The code only branches here once xm < tolbrent.
 2    continue

C     Final output
      if (verbose) then
         call output('write_fin', gamlist(its), det, fac, foundzero)
         call output('close', gamlist(its), det, fac, foundzero)
      end if

      end



C***********************************************************************


      subroutine output(action, gam, d1, fac, foundzero)

      implicit none

      character action*(*)
      double precision d1, gam, fac

      integer nymax, j
      include '../nymax.inc'
C     $   y(nymax), pip(nymax), psip(nymax), junk(nymax), xp
      logical CPexists, foundzero

      if (action.eq.'open') then
         inquire(file='pert_junk/status.junk', exist=CPexists)
         if (CPexists) then 
            open(unit=14, position='append', 
     $            file='pert_junk/doutdin.junk', status='old')
            open(unit=15, position='append',
     $            file='pert_junk/detofgam.junk', status='old')
            open(unit=16, position='append', 
     $            file='pert_junk/status.junk', status='old')

         else
            call system('rm -rf pert_junk; mkdir pert_junk')

            open(unit=14, file='pert_junk/doutdin.junk', status='new')
            open(unit=15, file='pert_junk/detofgam.junk', status='new')
            open(unit=16, file='pert_junk/status.junk', status='new')
         end if
      else if (action.eq.'write') then

         write(15,'(G24.16,G24.16,G24.16,L1)') gam, d1, fac, foundzero
         call flush(15)

      else if (action.eq.'write_fin') then

         open(unit=10, file='gam.out', status='unknown')
         write(10,99) gam
         close(10)

      else if (action.eq.'close') then
         do j=11,16
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
      write(16,*) its, ivar
      call flush(16)

 99   format(G24.16)

      return
      end

      
