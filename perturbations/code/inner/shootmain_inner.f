C***********************************************************************
C
C     Main program for calculating the Choptuik spacetime
C     (inner patch), modified version of code by CG & JM 
C
C     Structure of the program
C
C        1.- Variable declarations
C        2.- Read parameters
C        3.- Generate grid
C        4.- Input free data or previous checkpointing
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
      double precision in0(n3max), out0(n3max),
     $      lamb, fc(nymax), Up(nymax), psic(nymax), gam
 
C     x related variables
      integer nx, ileft, imid, iright, i, outevery, method,
     $     nxmax, itsreach, nleft, nright
      parameter(nxmax = 100001)
      double precision xc, xm, xp, xleft, xmid, xright,
     $     refine, crit, xxp(nxmax), dx

C     Background variables
      double precision Delta,  
     $     uB(nymax, nxmax), vB(nymax,nxmax), 
     $     fB(nymax,nxmax), ia2B(nymax, nxmax),
     $     fcB(nymax), psicB(nymax), UpB(nymax), in0B(n3max)

C     Newton method variables
      integer its, indx(n3max), ivar, maxits
      parameter(maxits=100)
      double precision in1(n3max), out1(n3max),
     $     doutdin(n3max,n3max), doutdinLU(n3max,n3max),
     $     eps, d, err, errmax, slowerr, fac,
     $     d1, gamlist(maxits), d1list(maxits),
     $     log_eig_sum, gamstep

C     Output
      integer itsread, ivarread
      logical verbose, debug, CPexists, foundzero

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
      read(10,*) refine
      read(10,*) method
      read(10,*) crit
      read(10,*) debug
      close(10)
     
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

C     Grid points:
C       ileft: leftmost point of the grid, corresponding to xleft
C       im: change of logarithmic axes, corresponding to xm
C       imid: junction point between left and right evolutions
C       iright: rightmost point of the grid, corresponding to xright

C     Point index
      ileft = 1
      iright = nleft + nright + 1

      nx = iright

C     Read in x-grid
      open(unit=10,file='grid.dat',status='old')
         do i=1,nx
            read(10,*) xxp(i), dx
         end do
      close(10)
      
C     Look for imid
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
     
C     *******************   
C     *** Background ****
C     *******************

C     Read in fields
      write(6,*) 'INFO: Reading in background'
      open(unit=10,file='uB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) uB(j,i)
            end do
         end do
      close(10)
      open(unit=10,file='vB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) vB(j,i)
            end do
         end do
      close(10)
      open(unit=10,file='fB.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) fB(j,i)
            end do
         end do
      close(10)
      open(unit=10,file='ia2B.dat',status='old')
         do i = ileft, iright
            do j = 1, ny
               read(10,*) ia2B(j,i)
            end do
         end do
      close(10)


C     Read in background initial data from out files. 
      open(unit=10, file='Delta.out', status='old')
        read(10,*) Delta
      close(10)
      open(unit=10, file='fcB.dat', status='old')
        do j=1,ny
          read(10,*) fcB(j)
        end do
      close(10)
      open(unit=10, file='psicB.dat', status='old')
        do j=1,ny
          read(10,*) psicB(j)
        end do
      close(10)    
      open(unit=10, file='UpB.dat', status='old')
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
      open(unit=10, file='delfc.dat', status='old')
      do j=1,ny
         read(10,*) fc(j)
      end do
      close(10)
      open(unit=10, file='delpsic.dat', status='old')
      do j=1,ny
         read(10,*) psic(j)
      end do
      close(10)    
      open(unit=10, file='delUp.dat', status='old')
      do j=1,ny
         read(10,*) Up(j)
      end do
      close(10)

C     Initialize output
      do its=1,maxits
         gamlist(its) = 0.d0
         d1list(its) = 0.d0
      end do
      
      gamlist(1) = gam
      gamstep = - 1.d-2
      

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
C        Read in doutdin.
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


C        Read in lamb and det
         open(unit=10, file='detofgam.junk', status='old')
         do its=1,itsread
               read(10,*) gamlist(its), d1list(its), fac, foundzero
         end do
         close(10)

      else

C     Brand new evolution.
         itsread = 0
         ivarread = 0
         fac = 1.d0
         foundzero = .FALSE.
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

C     Force log output
      call flush(6)

C     ******************
C     **** Shooting ****
C     ******************

      if (verbose)
     $  call output('open', gam, fac, d1, foundzero)
      
C     Shooting
      itsreach = 1
     

C     Iterations for finding zero of det(A)(lamb)
      do its=itsread,maxits-1
         
C        Set next value of lambda     
C        Compute new lamb from older two by secant method
C        or make another normal step if zero was not found yet
         if (its.gt.0) then                        
            if (foundzero) then
               gamlist(its+1) = (gamlist(its-1)*d1list(its) 
     $                           - gamlist(its)*d1list(its-1))
     $                        /(d1list(its) - d1list(its-1))            
            else
               gamlist(its+1) = gamlist(its) + gamstep
            end if
         end if
         
         lamb = 1.d0 / gamlist(its+1)
         
         write(6,*) 'gamma now = ', gamlist(its+1)

         
         write(6,*) 'INFO: Shooting iteration ', its
         call shootlin(ny, nx, d, lamb, ileft, iright, imid, 0,
     $      n3, in0, in0B, out0, outevery, xxp, method, 
     $      crit, debug, itsreach, 
     $      uB, vB, fB, ia2B)
         
         
C        Calculate mismatch (just for comparison, not important here)
         err = 0.d0
         do j=1,n3
            err = err + out0(j)**2
         end do
         err = sqrt(err/n3)
         write(*,*) 'Mismatch: ', err
         
C        Debug stop
         if (debug) stop
           
         
C        Compute Jacobian
         do ivar=ivarread+1,n3

C           Store input free data
            do j=1,n3
               in1(j) = in0(j)
            end do
C           Change variable ivar
            in1(ivar) = in1(ivar) + eps

C           Output info
            write(6,*) '  Max its reached:', itsreach, '. Varying', ivar
            call flush(6)

C           Shoot new free data
            itsreach = 1
            call shootlin(ny, nx, d, lamb, ileft, iright, imid, ivar,
     $           n3, in1, in0B, out1, 0, xxp, method, crit, debug,  
     $           itsreach, uB, vB, fB, ia2B)

C           Calculate derivative with respect to variable ivar
            do j=1,n3
               doutdin(j,ivar) = (out1(j) - out0(j)) / eps
            end do
            
            if (verbose) 
     $         call output_step(n3, its+1, ivar, doutdin(1,ivar))

         end do
      

C        For the fixed step iterations an exponent is determined which scales the 
C        matrix in order to render Det(doutdin*10^expon) of order 1
         if (foundzero) then
C           Compute determinant
C           LU decomposition of doutdin, stored in the same space
            call ludcmp(doutdin, n3, n3max, indx, d1)
            do j = 1, n3
               d1 = d1 * doutdin(j,j) * fac
            end do
         else
            call ludcmp(doutdin, n3, n3max, indx, d1)

            log_eig_sum = 0.d0

            do j = 1, n3
               log_eig_sum = log_eig_sum + log10(abs(doutdin(j,j)))
            end do

            fac = 10.d0**( - log_eig_sum / dble(n3))
            
            do j = 1, n3
               d1 = d1 * doutdin(j,j) * fac
            end do

            write(6,*) '(fac, det) = (', fac, d1, ' )'

         end if

         if ((its.gt.0).and.(.not.foundzero)) then
            if ((d1 / d1list(its)) .lt. 0.0d0) then
            
               foundzero = .TRUE.
               write(6,*) 'Found zero, switching to secant method.'
               
            end if
         end if

         d1list(its+1) = d1
       
         write(6,*) '(gamma, det) = (', gamlist(its+1), ' ,', d1,' )'

         
C        Write current step result to file
         if (verbose) then
               call output('write', 1.d0 / lamb, fac, d1, foundzero)
         end if
         
         ivarread = 0
         
C        If determinant is small enough, branch out
C         if (abs(d1list(its+1)-d1list(its)) .lt. 1.d-5) goto 1
C         if (its.eq.7) goto 1

      end do

C     The code only branches here once d1 < errmax.
 1    continue

C     Final output
      if (verbose) then
         call output('write_fin', 1.d0 / lamb, fac, d1, foundzero)
         call output('close', 1.d0 / lamb, fac, d1, foundzero)
      end if

      end



C***********************************************************************


      subroutine output(action, gam, fac, d1, foundzero)

      implicit none

      character action*(*)
      double precision fac, d1, gam

      integer nymax, j
      include '../nymax.inc'
C     $   y(nymax), pip(nymax), psip(nymax), junk(nymax), xp
      logical CPexists, foundzero

      if (action.eq.'open') then
         inquire(file='status.junk', exist=CPexists)
         if (CPexists) then 
            open(unit=14, position='append', 
     $            file='doutdin.junk', status='old')
            open(unit=15, file='detofgam.junk', 
     $            position='append', status='old')
            open(unit=16, position='append', 
     $            file='status.junk', status='old')

         else
            open(unit=14, file='doutdin.junk', status='new')
            open(unit=15, file='detofgam.junk', status='new')
            open(unit=16, file='status.junk', status='new')
         end if
      else if (action.eq.'write') then

         write(15,*) gam, d1, fac, foundzero
         call flush(15)

      else if (action.eq.'write_fin') then

C         open(unit=10, file='lamb.out', status='unknown')
C         write(10,99) lamb
C         close(10)

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

      
