C***********************************************************************
C
C    Subroutine of shootmain_inner
C         
C
C    Given a tolerance tol, free data is shot to the point where psi
C    takes half the value of its maximum in the center. 
C    Adaptive stepsize is used to generate the grid in x 
C
C    OUTPUT: xxp, nleft, nright, xmid
C
C***********************************************************************

      subroutine gen_grid(ny, d, xleft, xright, xmid,
     $     n3, in, out, prec_irk, tol, xxp, nleft, nright, 
     $     dxini, debug)

      implicit none
      
      integer ny, nleft, nright, n3, nxmax, maxits, nymax
      include '../nymax.inc'
      parameter(nxmax = 30001, maxits = 20)

      double precision in(n3), out(n3), xxp(nxmax), prec_irk, d, dxin, 
     $                 dxguess, tol, xin, xout, 
     $                 xleft, xright, xmid, dxini, psicmax, psinowmax

      double precision fc(nymax), psic(nymax), Up(nymax), Delta,
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     junk(nymax),
     $     ypleft(nymax), ypright(nymax), mismatch(nymax),
     $     xxp_temp(nxmax) 
      integer i, j
      logical debug

      external derivs

C     Safety checks.
      if (4*n3 .ne. 3*ny) stop 'n3 .ne. 3/4 ny'

C     Unpack free data.
      Delta = in(ny/2+3)
      in(ny/2+3) = 0.d0
      call myunpack(n3, in, ny, Up, psic, fc)
      in(ny/2+3) = Delta

C     Determine maximum value of psic for fixing xmid at 
C     half the value of that maximum
      psicmax = 0.d0
      do i=1,ny
            psicmax = max(psicmax,psic(i))
      end do

      write(6,*) 'Psicmax = ', psicmax

C     *************************
C     Shoot from ileft to imid.
C     *************************
      call leftdata(ny, d, xleft, Delta, fc, psic, ypleft, 
     $              debug, .true.)

      dxin = dxini
      xin = xleft

      write(6,*) 'Generating left grid...'
      do i = 1, nxmax
         
C        write(6,*) i, xin
        xxp(i) = xin
        call irk2_stepper(ny, d, xin, ypleft, xout, 
     $         ypleft, prec_irk, tol, dxin, dxguess, maxits, derivs)
    
        dxin = dxguess
        xin = xout

C       Compute psimax at new step and check whether it is already below
C       half maximum -> determines xmid
        call fieldsfromy(ny, d, ypleft, xout,
     $     u, v, f, ia2, Delta,
     $     junk, junk, junk)

        psinowmax = 0.d0
        do j=1,ny
            psinowmax = max(psinowmax,(v(j)-u(j))/(2.d0*xout**2))
        end do

        if (psinowmax.lt.psicmax/2.d0) then
            xmid = xout
            xxp(i+1) = xmid
            nleft = i
            exit
        end if

C        version for fixed predefined xmid
C        if (xout.gt.xmid) then
C            xxp(i+1) = xmid
C            nleft = i
C            exit
C        end if
         
      end do
      
C     **************************      
C     Shoot from iright to imid.
C     ************************** 
      write(6,*) 'Generating right grid...'
      call rightdata(ny, d, xright, Delta, Up, ypright, 
     $       debug, .true.)

      dxin = - dxini
      xin = xright

      do i = 1, nxmax

C        write(6,*) i, xin
        xxp_temp(i) = xin
        call irk2_stepper(ny, d, xin, ypright, xout, 
     $         ypright, prec_irk, tol, dxin, dxguess, maxits, derivs)
      
        dxin = dxguess
        xin = xout

        if (xout.lt.xmid) then
            nright = i
            exit
        end if
      end do

C      write(6,*) nleft, nright
      do i=1,nright
      xxp(1+nleft+i) = xxp_temp(nright - i + 1)
      end do


C     Calculate mismatch. IRRELEVANT HERE
      do j=1,ny
         mismatch(j) = ypright(j) - ypleft(j)
      end do
      if (mismatch(5).ne.0.d0) stop 'Error in Delta'
      mismatch(5) = Delta

      call fieldsfromy(ny, d, mismatch, xmid,
     $     u, v, f, ia2, Delta,
     $     junk, junk, junk)
     
      call mypack(n3, out, ny, u, v, f)

      return
      end
