C***********************************************************************
C
C    Subroutine of shootmain_inner
C                  modified version of code by CG & JM 
C
C    Shoots a given free data set from xleft to xmid, and from xright
C    to xmid and computes mismatch.

C
C***********************************************************************

      subroutine shoot_to_file(ny, nx, d, ileft, iright, imid,
     $     n3, in, out, xxp, tstep, prec_irk, debug, itsreach)

      implicit none
      integer ny, nx, ileft, iright, imid, n3, itsreach
      double precision in(n3), out(n3), xxp(nx), prec_irk, d
      logical debug

      integer nymax, nxmax, maxits
      include '../nymax.inc'
      parameter(nxmax = 300001, maxits = 100)

      double precision fc(nymax), psic(nymax), Up(nymax), Delta,
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     junkxi(nymax),
     $     yp(nymax, nxmax),
     $     mismatch(nymax), ytest(nymax)
      integer i, j, ifile, tstep
      logical repeat
      external derivs

      save yp


C     Safety checks.
      if (nx .gt. nxmax) stop 'nx gt nxmax in shoot1'
      if (4*n3 .ne. 3*ny) stop 'n3 .ne. 3/4 ny'

C     Unpack free data.
      Delta = in(ny/2+3)
      in(ny/2+3) = 0.d0
      call myunpack(n3, in, ny, Up, psic, fc)
      in(ny/2+3) = Delta
             
             
C     Shoot from ileft to imid.
      call leftdata(ny, d, xxp(ileft), Delta, fc, psic, 
     $              yp(1,ileft), debug, tstep)
      
      do i = ileft, imid-1
            if (MOD(i,500).eq.0) then
               write(6,*) i
            end if
            call implicit_step(ny, d, xxp(i), yp(1,i), xxp(i+1),  
     $        yp(1,i+1), derivs, 2, prec_irk, maxits, itsreach,
     $        repeat)
      
      end do
      
C     Shoot from iright to imid+1. 
      call rightdata(ny, d, xxp(iright), Delta, Up, yp(1,iright),
     $  debug, tstep)
      do i = iright, imid+2, -1
       
            call implicit_step(ny, d, xxp(i), yp(1,i), xxp(i-1),         
     $         yp(1,i-1), derivs, 2, prec_irk, maxits, itsreach,
     $         repeat)
    
      end do


C     Shoot one further from right to imid and compare (for computing mismatch)
      call implicit_step(ny, d, xxp(imid+1), yp(1,imid+1), xxp(imid),
     $    ytest, derivs, 2, prec_irk, maxits, itsreach, repeat)
      

      do j=1,ny
         mismatch(j) = ytest(j) - yp(j,imid)
      end do
      
      if (mismatch(5).ne.0.d0) stop 'Error in Delta'
      mismatch(5) = Delta
     

      call fieldsfromy(ny, d, mismatch, xxp(imid),
     $     u, v, f, ia2, Delta,
     $     junkxi, junkxi, junkxi)
     
     
      call mypack(n3, out, ny, u, v, f)
      
C     Output
      open(unit=80,file='bg_data/uB.dat',status='new')
      open(unit=81,file='bg_data/vB.dat',status='new')
      open(unit=82,file='bg_data/fB.dat',status='new')
      open(unit=83,file='bg_data/ia2B.dat',status='new')
         
      write(6,*) 'INFO: Writing to file...'

         do i=1,nx       
            call fieldsfromy(ny, d, yp(1,i), xxp(i), 
     $           u, v, f, ia2, Delta,
     $           junkxi, junkxi, junkxi)
            if (tstep.eq.1) then
               do j=1,ny
                  write(80,99) u(j)
                  write(81,99) v(j)
                  write(82,99) f(j)
                  write(83,99) ia2(j)
               end do
            else
               call changeresolution(ny, ny/2, u, u)
               call changeresolution(ny, ny/2, v, v)
               call changeresolution(ny, ny/2, f, f)
               call changeresolution(ny, ny/2, ia2, ia2)
               do j=1,ny/2
                  write(80,99) u(j)
                  write(81,99) v(j)
                  write(82,99) f(j)
                  write(83,99) ia2(j)
               end do
            end if
            do ifile=80,83
               write(ifile,*) ' '
            end do
         
         end do
         
         do ifile=80,83
            close(ifile)
         end do

99    format(G24.16)

      return
      end
