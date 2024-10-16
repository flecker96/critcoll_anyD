C***********************************************************************
C
C    Subroutine of shootmain_inner
C                  modified version of code by CG & JM 
C
C    Shoots a given free data set from xleft to xmid, and from xright
C    to xmid and computes mismatch.
C
C
C***********************************************************************

      subroutine shoot1(ny, nx, d, ileft, iright, imid, ivar,
     $     n3, in, out, outevery, xxp, prec_irk, debug, 
     $     itsreach)

      implicit none
      integer ny, nx, ileft, iright, imid, ivar, n3, itsreach
      double precision in(n3), out(n3), xxp(nx), prec_irk, d
      logical debug

      integer nymax, nxmax, maxits
      include '../nymax.inc'
      parameter(nxmax = 30001, maxits = 100)

      double precision fc(nymax), psic(nymax), Up(nymax), Delta,
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     junk(nymax),
     $     ypleft(nymax), ypright(nymax), mismatch(nymax),
     $     ypleftbase(nymax), yprightbase(nymax)
      integer i, j, outevery
      logical evl, evr, repeat, printtayl

      external derivs, implicit_step

      save ypleftbase, yprightbase

C     Safety checks.
      if (nx .gt. nxmax) stop 'nx gt nxmax in shoot1'
      if (4*n3 .ne. 3*ny) stop 'n3 .ne. 3/4 ny'

C     Unpack free data.
      Delta = in(ny/2+3)
      in(ny/2+3) = 0.d0
      call myunpack(n3, in, ny, Up, psic, fc)
      in(ny/2+3) = Delta
      
             
C     Output
      if (outevery.ne.0) then
         open(unit=79,file='xaxis.junk',status='unknown')
         open(unit=80,file='u.junk',status='unknown')
         open(unit=81,file='v.junk',status='unknown')
         open(unit=82,file='f.junk',status='unknown')
         open(unit=83,file='ia2.junk',status='unknown') 
      end if

C     Variables junkc, junkf1 and junkf2 are irrelevant in the inner
C     patch.

      if (ivar.eq.0) then
         evl = .TRUE.
         evr = .TRUE.
         printtayl = .true.
      else if (ivar.le.n3/3) then
         evl = .FALSE.
         evr = .TRUE.
         printtayl = .false.
      else if (ivar.eq.2*n3/3+3) then
         evl = .TRUE.
         evr = .TRUE.
         printtayl = .false.
      else
         evl = .TRUE.
         evr = .FALSE.
         printtayl = .false.
      end if

      if (evl) then
C     Shoot from ileft to imid.
      call leftdata(ny, d, xxp(ileft), Delta, fc, psic, ypleft, 
     $               debug, printtayl)
      do i = ileft, imid-1
C     Two nested ifs because some compilers check both questions,
C     giving an error when outevery=0.
         if (outevery.ne.0) then
C         if (i.eq.imid) then
         if (mod(i,outevery).eq.0) then
            write(6,*) 'writing to file', i, xxp(i)
            call fieldsfromy(ny, d, ypleft, xxp(i),
     $           u, v, f, ia2, Delta,
     $           junk, junk, junk)
            write(79,99) xxp(i)
            do j=1,ny
               write(80,99) u(j)
               write(81,99) v(j)
               write(82,99) f(j)
               write(83,99) ia2(j)
            end do
            
         end if
         end if
         call implicit_step(ny, d, xxp(i), ypleft, xxp(i+1),  
     $        ypleft, derivs, 2, prec_irk, maxits, itsreach, 
     $        repeat)
         if (repeat) then 
           write(6,*) 'implicit step did not converge at x = ', xxp(i)
           stop
         end if
      end do
      else
         do j=1,ny
            ypleft(j) = ypleftbase(j)
         end do
      end if

      if (evr) then
C     Shoot from iright to imid. 
      call rightdata(ny, d, xxp(iright), Delta, Up, ypright, 
     $    debug, printtayl)
      do i = iright, imid+1, -1
C     Two nested ifs because some compilers check both questions,
C     giving an error when outevery=0.
C         write(6,*) xxp(i)
         if (outevery.ne.0) then
         if (mod(i,outevery).eq.0) then
C         if (i.eq.imid) then
            call fieldsfromy(ny, d, ypright, xxp(i),
     $           u, v, f, ia2, Delta,
     $           junk, junk, junk)
            write(79,99) xxp(i)
            do j=1,ny
               write(80,99) u(j)
               write(81,99) v(j)
               write(82,99) f(j)
               write(83,99) ia2(j)
            end do

         end if
         end if
         call implicit_step(ny, d, xxp(i), ypright, xxp(i-1),
     $        ypright, derivs, 2, prec_irk, maxits, itsreach, 
     $        repeat)
         if (repeat) then 
           write(6,*) 'implicit step did not converge at x = ', xxp(i)
           stop
         end if
      end do
      else
         do j=1,ny
            ypright(j) = yprightbase(j)
         end do
      end if

C     Save base results
      if (ivar.eq.0) then
         do j=1,ny
            ypleftbase(j) = ypleft(j)
            yprightbase(j) = ypright(j)
         end do
      end if
      
C     Calculate mismatch.
      do j=1,ny
         mismatch(j) = ypright(j) - ypleft(j)
      end do
      if (mismatch(5).ne.0.d0) stop 'Error in Delta'
      mismatch(5) = Delta

      call fieldsfromy(ny, d, mismatch, xxp(imid),
     $     u, v, f, ia2, Delta,
     $     junk, junk, junk)
     
      
      call mypack(n3, out, ny, u, v, f)

      do j=79,83
         close(j)
      end do

99    format(F24.16)

      return
      end
