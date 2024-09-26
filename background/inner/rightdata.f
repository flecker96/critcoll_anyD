      subroutine rightdata(ny, d, x, Delta, Up, y, debug)

      implicit none
      integer ny
      double precision x, d, Delta, Up(ny), y(ny)
      logical debug

      integer nymax, j
      include '../nymax.inc'
      double precision xp
      double precision coeff1(nymax), coeff2(nymax),
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     pi(nymax), psi(nymax),
     $     u0(nymax), v0(nymax), f0(nymax), ia20(nymax),
     $     u1(nymax), v1(nymax), f1(nymax), ia21(nymax),
     $     u2(nymax), v2(nymax), f2(nymax), ia22(nymax),
     $     du0dtau(nymax), d2u0dtau2(nymax), junkxi(nymax)

      if (debug) write(6,*) 'rightdata debug at x=', x

      xp = 1.d0

C *** Order 0 ***

      do j=1,ny
         f0(j) = xp
         u0(j) = Up(j)
      end do
      call fourdiff1(ny, u0, du0dtau, Delta)
      call fourdiff1(ny, du0dtau, d2u0dtau2, Delta)
c     call fourdiff1smoothe(ny, U0, dU0dtau, Delta)
c     call fourdiff1smoothe(ny, dU0dtau, d2U0dtau2, Delta)
      
    
      
C     Solve a linear inhomogeneous ODE for ia20 
      do j=1,ny
         coeff1(j) = 3.d0 - d - (d - 2.d0)**3 * u0(j)**2 / 4.d0 
         coeff2(j) = d - 3.d0
      end do
      call inhom(ia20, coeff1, coeff2, Delta, ny)

C     Solve a linear inhomogeneous ODE for v0. 
      do j=1,ny
         coeff1(j) = (6.d0 - 2.d0*d + (d - 2.d0) * ia20(j)) 
     $                / (2.d0*ia20(j))
         coeff2(j) = u0(j) * (d - 2.d0) / 2.d0
      end do
      call inhom(v0, coeff1, coeff2, Delta, ny)

C *** Order 1 ***

C     f1, u1, ia21 are given algebraically
      do j=1,ny
         f1(j) = 3.d0 - d + (d - 3.d0) / ia20(j)   
         ia21(j) = - ia20(j) * (8.d0*(d - 3.d0) 
     $                          + (d - 2.d0)**3 * (u0(j)**2 + v0(j)**2))
     $             / 8.d0 + d - 3.d0
         u1(j) = (d - 2.d0 - 2.d0 * (d - 3.d0) / ia20(j)) * u0(j) / 4.d0
     $            + (d - 2.d0) * v0(j) / 4.d0 - du0dtau(j) / 2.d0   
      end do

C     Solve a linear inhomogeneous ODE for v1. 
      do j=1,ny
         coeff1(j) = (6.d0 - 2.d0*d + d * ia20(j) - 2.d0*f1(j)*ia20(j))
     $                / (2.d0*ia20(j)) 
         coeff2(j) = (2.d0 * (3.d0 - d)*(f1(j) - 1.d0) * ia20(j) * v0(j)
     $                + 2.d0 * (d - 3.d0) * ia21(j) * v0(j)
     $                + (d - 2.d0)*ia20(j)**2 * ((f1(j) - 1.d0) * u0(j)
     $                   + u1(j) + (f1(j) - 1.d0) * v0(j))) 
     $               / (2.d0*ia20(j)**2)    
      end do
      call inhom(v1, coeff1, coeff2, Delta, ny)

C *** Order 2 ***

C     f2, u2, ia22 are given algebraically
      do j=1,ny
         f2(j) =  (d - 3.d0) * ((f1(j) - 1.d0)*(1.d0 - ia20(j))*ia20(j)
     $                          - ia21(j)) / (2.d0 * ia20(j)**2)    
         ia22(j) = (3.d0 - d - (ia21(j) - ia20(j))*(8.d0*(d - 3.d0)
     $                 + (d - 2.d0)**3 * (u0(j)**2 + v0(j)**2)) / 8.d0
     $              - (d - 2.d0)**3 * ia20(j)
     $                * (u0(j)*u1(j) + v0(j)*v1(j) ) / 4.d0 ) / 2.d0
         u2(j) = ((4.d0*(3.d0 - d)*(d - 6.d0 + f1(j))*ia20(j)  
     $             + (d - 2.d0)*(d - 8.d0 + 2.d0*f1(j))*ia20(j)**2
     $             + 4.d0*(d - 3.d0)*(d - 3.d0 + 2.d0*ia21(j))) * u0(j)
     $           - (d - 3.d0)*(d - 2.d0)**3 * ia20(j) * u0(j)**3
     $           + (4.d0*(6.d0 - 2.d0*d + (d - 2.d0)*ia20(j))*u1(j)
     $              + (d - 2.d0) * v0(j) * (6.d0 - 2.d0*d  
     $                               + (d - 8.d0 + 2.d0*f1(j))*ia20(j))
     $              - 8.d0*ia20(j)*v1(j) + 4.d0*d*ia20(j)*v1(j)
     $              - 12.d0*du0dtau(j) + 4.d0*d*du0dtau(j)
     $              + 8.d0*ia20(j)*du0dtau(j)- 2.d0*d*ia20(j)*du0dtau(j)
     $              + 4.d0*f1(j)*ia20(j)*du0dtau(j)
     $              + 4.d0*ia20(j)*d2u0dtau2(j)) * ia20(j) ) 
     $           / (32.d0 * ia20(j)**2)
      end do
      
C     Solve a linear inhomogeneous ODE for v2. 
      do j=1,ny
         coeff1(j) = 1.d0 + d / 2.d0 - 2.d0*f1(j) + (3.d0 - d) / ia20(j)
         
         coeff2(j) = (2.d0*(3.d0 - d)*ia21(j)**2 *v0(j)+ 2.d0*(d - 3.d0)
     $                * ia20(j)**2 * ((f1(j) - 1.d0 - f2(j)) * v0(j)
     $                                - (f1(j) - 1.d0)*v1(j))
     $                + ((d - 2.d0)*((1.d0 - f1(j) + f2(j)) * u0(j)
     $                     + (f1(j) - 1.d0) * u1(j) + u2(j)
     $                     + (1.d0 - f1(j) + f2(j)) * v0(j))
     $                  + ((d - 2.d0)*(f1(j) - 1.d0) - 2.d0*f2(j))*v1(j)
     $                  ) * ia20(j)**3
     $                  + 2.d0*(d - 3.d0)*ia20(j)*(ia22(j)*v0(j)  
     $                  + ia21(j) * ((f1(j) - 1.d0)*v0(j) + v1(j) ) )
     $                   ) / (2.d0 * ia20(j)**3)
      end do
      call inhom(v2, coeff1, coeff2, Delta, ny)

C     Put together the expansion around x=xp=1.
      do j=1,ny
         f(j)  = f0(j) + (x-xp) * f1(j) + (x-xp)**2 * f2(j)
         ia2(j) = ia20(j) + (x-xp) * ia21(j) + (x-xp)**2 * ia22(j)
         u(j)  = u0(j)+ (x-xp) * u1(j) + (x-xp)**2 * u2(j)
         v(j)  = v0(j)+ (x-xp) * v1(j) + (x-xp)**2 * v2(j)
C         psi(j)  = ( V(j) - U(j) ) / 2.d0 / x**2
C         pi(j) = ( V(j) + U(j) ) / 2.d0 / x
      end do
     
C               open(unit=10, file='transfer.junk', status='new')
C             open(unit=11, file='transfer2.junk', status='new')
C             open(unit=12, file='transfer3.junk', status='new')
C               do j=1,ny
C                  write(10,*) u(j), f(j), v(j)
C                write(11,*) psil(i), psir(i)
C                  write(12,*) fl(i), fr(i)
C               end do
C               close(10) 
C              close(11)
C               close(12)
C              stop

     
      do j=1,ny
         junkxi(j) = 0.d0
      end do

      call yfromfields(ny, u, v, f, y)

      y(5) = Delta

C     Write out taylor coefficients      
      if (debug) then
         
            open(unit=11,file='R_taylor/fright.junk',status='new')
            open(unit=12,file='R_taylor/uright.junk',status='new')
            open(unit=13,file='R_taylor/vright.junk',status='new')
            open(unit=14,file='R_taylor/ia2right.junk',status='new')
      
               do j=1,ny
                  write(11,*) f0(j), f1(j), f2(j)
                  write(12,*) u0(j), u1(j), u2(j)
                  write(13,*) v0(j), v1(j), v2(j)
                  write(14,*) ia20(j), ia21(j), ia22(j)
               end do
      
            do j=11,14
               close(j)
            end do
      end if

99    format(F24.16)

      return
      end

