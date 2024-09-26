C***********************************************************************
C
C    Subroutine of shootmain_inner
C                  modified version of code by CG & JM 
C
C    Shoots a given free data set from xleft to xmid, and from xright
C    to xmid and computes mismatch.

C
C***********************************************************************

      subroutine shootlin(ny, nx, d, lamb, ileft, iright, imid, ivar,
     $     n3, in, inB, out, outevery, xxp, method, crit, debug, 
     $     itsreach, uB, vB, fB, ia2B)

      implicit none
      integer ny, nx, ileft, iright, imid, ivar, n3, method, itsreach
      double precision in(n3), out(n3), xxp(nx), crit, d, lamb
      logical debug

      integer nymax, nxmax, maxits
      include '../nymax.inc'
      parameter(nxmax = 100001, maxits = 100)

      double precision fc(nymax), psic(nymax), Up(nymax), Delta,
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     ypleft(nymax), ypright(nymax), mismatch(ny),
     $     ypleftbase(nymax), yprightbase(nymax), delfc(nymax),
     $     delpsic(nymax), delUp(nymax), 
     $     uB(nymax,nxmax), vB(nymax,nxmax),  inB(ny),
     $     fB(nymax,nxmax), ia2B(nymax,nxmax), c(2),
     $     du(nymax), dv(nymax), df(nymax), dia2(nymax),
     $     uBcoll(nymax,2), vBcoll(nymax,2), fBcoll(nymax,2),
     $     ia2Bcoll(nymax,2), junk(nymax)
      integer i, j, k, outevery
      logical evl, evr
      external derivslin

      save ypleftbase, yprightbase

C     Collocation points
      c(1) = 0.5d0 - 0.5d0 / sqrt(3.d0)
      c(2) = 0.5d0 + 0.5d0 / sqrt(3.d0)

C     Safety checks.
      if (nx .gt. nxmax) stop 'nx gt nxmax in shoot1'
      if (4*n3 .ne. 3*ny) stop 'n3 .ne. 3/4 ny'

C     Unpack free data.
      call myunpack(n3, in, ny, delUp, delpsic, delfc)

C     Unpack background data.
      Delta = inB(ny/2+3)
      inB(ny/2+3) = 0.d0
      call myunpack(n3, inB, ny, Up, psic, fc)
      inB(ny/2+3) = Delta
      
             
C     Output
      if (outevery.ne.0) then
         open(unit=79,file='xaxis.junk',status='unknown')
         open(unit=80,file='u.junk',status='unknown')
         open(unit=81,file='v.junk',status='unknown')
         open(unit=82,file='f.junk',status='unknown')
         open(unit=83,file='ia2.junk',status='unknown') 
      end if

C     Shoot only from the side which is varied
      if (ivar.eq.0) then
         evl = .TRUE.
         evr = .TRUE.
      else if (ivar.le.n3/3) then
         evl = .FALSE.
         evr = .TRUE.
      else
         evl = .TRUE.
         evr = .FALSE.
      end if
      
C      evl = .FALSE.
      if (evl) then
C     Shoot from ileft to imid.
      call leftdatalin(ny, d, xxp(ileft), Delta, fc, psic,
     $              delfc, delpsic, ypleft, lamb, debug)  
C      call leftdatalin(ny, d, xxp(imid), Delta, fc, psic,
C     $              delfc, delpsic, ypleft, lamb, debug)  

      do i = ileft, imid-1
C      do i = imid, ileft+1, -1
C     Two nested ifs because some compilers check both questions,
C     giving an error when outevery=0.
         if (outevery.ne.0) then
C         if (i.eq.imid-1) then
         if ((mod(i,outevery).eq.0).or.(i.eq.ileft)) then
C            call leftdatalin(ny, d, xxp(i), Delta, fc, psic,
C     $              delfc, delpsic, ypleft, lamb, debug)  
       
            call fieldsfromylin(ny, d, ypleft, xxp(i), 
     $             u, v, f, ia2, Delta, lamb,
     $             junk, junk, junk,
     $             uB(1,i), vB(1,i), fB(1,i), ia2B(1,i))

            write(79,99) xxp(i)
            do j=1,ny
               write(80,99) u(j)
               write(81,99) v(j)
               write(82,99) f(j)
               write(83,99) ia2(j)
            end do
                   
         end if
         end if
         
C           Compute background at collocation points from uB(1,i), uB(1,i+1) etc.
C           Simple start: linear interpolation            
            do k=1,ny
                  du(k) = uB(k,i+1) - uB(k,i)
                  dv(k) = vB(k,i+1) - vB(k,i)
                  df(k) = fB(k,i+1) - fB(k,i)
                  dia2(k) = ia2B(k,i+1) - ia2B(k,i)
C                  du(k) = uB(k,i-1) - uB(k,i)
C                  dv(k) = vB(k,i-1) - vB(k,i)
C                  df(k) = fB(k,i-1) - fB(k,i)
C                  dia2(k) = ia2B(k,i-1) - ia2B(k,i)
               do j=1,2   
                  uBcoll(k,j) = du(k) * c(j) + uB(k,i)
                  vBcoll(k,j) = dv(k) * c(j) + vB(k,i)
                  fBcoll(k,j) = df(k) * c(j) + fB(k,i)
                  ia2Bcoll(k,j) = dia2(k) * c(j) + ia2B(k,i)
               end do
            end do
            
            call implicit_steplin(ny, d, Delta, lamb, 
     $        xxp(i), ypleft, xxp(i+1), ypleft, 
     $        method, crit, maxits, itsreach,
     $        uBcoll, vBcoll, fBcoll, ia2Bcoll)
      
      end do
      else
         do j=1,ny
            ypleft(j) = ypleftbase(j)
         end do
      end if
        
      

      if (evr) then
C     Shoot from iright to imid. 
      call rightdatalin(ny, d, xxp(iright), Delta, Up, delUp,
     $                  ypright, lamb, debug)
      do i = iright, imid+1, -1
C     Two nested ifs because some compilers check both questions,
C     giving an error when outevery=0.
      if (outevery.ne.0) then
          if ((mod(i,outevery).eq.0).or.(i.eq.iright))  then
C             call rightdatalin(ny, d, xxp(i), Delta, Up, delUp,
C     $                  ypright, lamb, debug)
             call fieldsfromylin(ny, d, ypright, xxp(i), 
     $             u, v, f, ia2, Delta, lamb,
     $             junk, junk, junk,
     $             uB(1,i), vB(1,i), fB(1,i), ia2B(1,i))
            write(79,99) xxp(i)
            do j=1,ny
               write(80,99) u(j)
               write(81,99) v(j)
               write(82,99) f(j)
               write(83,99) ia2(j)
            end do         
         end if
         end if
         
C           Compute background at collocation points from uB(1,i), uB(1,i+1) etc.
C           Simple start: linear interpolation  
            do k=1,ny
                  du(k) = uB(k,i-1) - uB(k,i)
                  dv(k) = vB(k,i-1) - vB(k,i)
                  df(k) = fB(k,i-1) - fB(k,i)
                  dia2(k) = ia2B(k,i-1) - ia2B(k,i)
                  
               do j=1,2   
                  uBcoll(k,j) = du(k) * c(j) + uB(k,i)
                  vBcoll(k,j) = dv(k) * c(j) + vB(k,i)
                  fBcoll(k,j) = df(k) * c(j) + fB(k,i)
                  ia2Bcoll(k,j) = dia2(k) * c(j) + ia2B(k,i)
               end do          
            end do 
         
        
            call implicit_steplin(ny, d, Delta, lamb, 
     $         xxp(i), ypright, xxp(i-1), ypright,        
     $         method, crit, maxits, itsreach,
     $         uBcoll, vBcoll, fBcoll, ia2Bcoll)
      
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

      call fieldsfromylin(ny, d, mismatch, xxp(imid), 
     $             u, v, f, ia2, Delta, lamb,
     $             junk, junk, junk,
     $             uB(1,imid), vB(1,imid), fB(1,imid), ia2B(1,imid))

      call mypack(n3, out, ny, u, v, f)

      do j=79,83
         close(j)
      end do

99    format(F24.16)

      return
      end
