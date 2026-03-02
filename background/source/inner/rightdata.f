      subroutine rightdata(ny, d, x, Delta, Up, y, debug, printtayl)

      implicit none
      integer ny
      double precision x, d, Delta, Up(ny), y(ny)
      logical debug, printtayl

      integer nymax, j
      include '../nymax.inc'
      double precision xp
      double precision coeff1(nymax), coeff2(nymax),
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     u0(nymax), v0(nymax), f0(nymax), ia20(nymax),
     $     u1(nymax), v1(nymax), f1(nymax), ia21(nymax),
     $     u2(nymax), v2(nymax), f2(nymax), ia22(nymax),
     $     u3(nymax), v3(nymax), f3(nymax), ia23(nymax),
     $     du0dtau(nymax), d2u0dtau2(nymax), 
     $     du2dtau(nymax)

C     Taylor variables     
      double precision f10max, f21max, u10max, u21max,
     $     v10max, v21max, 
     $     f0norm, f1norm, f2norm, u0norm, u1norm, u2norm,
     $     v0norm, v1norm, v2norm

      xp = 1.d0

C *** Order 0 ***

      do j=1,ny
         f0(j) = xp
         u0(j) = Up(j)
      end do
      call fourdiff1(ny, u0, du0dtau, Delta)
      call fourdiff1(ny, du0dtau, d2u0dtau2, Delta)
      
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
      
      call fourdiff1(ny, u2, du2dtau, Delta)

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


C *** Order 3 ***

      do j=1,ny
            f3(j) = ((-3.d0 + d)*((1.d0 - f1(j) + f2(j))*ia20(j)**2 + 
     $             (-1.d0 + f1(j) - f2(j))*ia20(j)**3 + ia21(j)**2 - 
     $             ia20(j)*((-1.d0 + f1(j))*ia21(j) + ia22(j))))/
     $         (3.d0*ia20(j)**3)
            ia23(j) = (-3.d0 + d + (-((ia20(j) - ia21(j) + ia22(j))*
     $                 (8.d0*(-3.d0 + d) + 
     $                   (-2.d0 + d)**3*(u0(j)**2 + v0(j)**2))) + 
     $              2.d0*(-2.d0 + d)**3*(ia20(j) - ia21(j))*
     $               (u0(j)*u1(j) + v0(j)*v1(j)) - 
     $              (-2.d0 + d)**3*ia20(j)*
     $               (u1(j)**2 + 2.d0*u0(j)*u2(j) + v1(j)**2 + 
     $                 2.d0*v0(j)*v2(j)))/8.d0)/3.d0
            u3(j) = -(16.d0*(-3.d0 + d)*ia21(j)**2*u0(j) + 
     $            4.d0*(-3.d0 + d)*ia20(j)*
     $             ((-3.d0 + d + f1(j)*(-3.d0 + d - 2.d0*ia21(j)) + 
     $                  6.d0*ia21(j) - 4.d0*ia22(j))*u0(j) - 
     $               4.d0*ia21(j)*u1(j)) - 
     $            (-3.d0 + d)*ia20(j)**2*
     $             (4.d0*(-10.d0 + d + f1(j)*(-1.d0 + d + f1(j)) - 
     $             2.d0*f2(j))*
     $                u0(j) + (-2.d0 + d)**3*(1.d0 + f1(j))*u0(j)**3 + 
     $               2.d0*(-4.d0*(-3.d0 + f1(j))*u1(j) - 8.d0*u2(j) + 
     $                  (1.d0 + f1(j))*
     $                   (-2.d0*du0dtau(j) + (-2.d0 + d)*v0(j)))) + 
     $            ia20(j)**3*(8.d0*du0dtau(j) - 2.d0*d*du0dtau(j) + 
     $               12.d0*du0dtau(j)*f1(j) - 2.d0*d*du0dtau(j)*f1(j) + 
     $               4.d0*du0dtau(j)*f1(j)**2 + 
     $               4.d0*d2u0dtau2(j)*(1.d0 + f1(j)) - 
     $               8.d0*du0dtau(j)*f2(j) + 
     $               (-2.d0 + d)*(-16.d0 + d + 
     $                  f1(j)*(2.d0 + d + 2.d0*f1(j)) - 
     $                  4.d0*f2(j))*u0(j) - 
     $               4.d0*(-2.d0 + d)*(-3.d0 + f1(j))*u1(j) 
     $                  + 16.d0*u2(j) - 
     $               8.d0*d*u2(j) + 32.d0*v0(j) - 18.d0*d*v0(j) + 
     $               d**2*v0(j) - 4.d0*f1(j)*v0(j) + 
     $               d**2*f1(j)*v0(j) - 4.d0*f1(j)**2*v0(j) + 
     $               2.d0*d*f1(j)**2*v0(j) + 8.d0*f2(j)*v0(j) - 
     $               4.d0*d*f2(j)*v0(j) - 24.d0*v1(j) + 12.d0*d*v1(j) + 
     $               8.d0*f1(j)*v1(j) - 4.d0*d*f1(j)*v1(j) 
     $                  + 16.d0*v2(j) - 
     $               8.d0*d*v2(j) + 
     $               16.d0*du2dtau(j)))/
     $          (96.d0*ia20(j)**3)
      end do

C     Solve a linear inhomogeneous ODE for v3. 
      do j=1,ny
         coeff1(j) = 2.d0 + d / 2.d0 - 3.d0*f1(j) + (3.d0 - d) / ia20(j)
         
         coeff2(j) = ((-3.d0 + d)*ia21(j)**3*v0(j))/ia20(j)**4 - 
     $      ((-3.d0 + d)*ia21(j)*
     $         (2.d0*ia22(j)*v0(j) + 
     $           ia21(j)*((-1.d0 + f1(j))*v0(j) + v1(j))))/
     $       ia20(j)**3 - ((-3.d0 + d)*
     $         ((-1.d0 + f1(j) - f2(j) + f3(j))*v0(j) + 
     $           (1.d0 - f1(j) + f2(j))*v1(j) + (-1.d0 + f1(j))*v2(j))
     $         )/ia20(j) + ((-2.d0 + d)*
     $          (-1.d0 + f1(j) - f2(j) + f3(j))*u0(j) - 
     $         (-2.d0 + d)*(-1.d0 + f1(j) - f2(j))*u1(j) + 2.d0*u2(j) - 
     $         d*u2(j) - 2.d0*f1(j)*u2(j) + d*f1(j)*u2(j) - 
     $         2.d0*u3(j) + d*u3(j) + 2.d0*v0(j) - d*v0(j) - 
     $         2.d0*f1(j)*v0(j) + d*f1(j)*v0(j) + 2.d0*f2(j)*v0(j) - 
     $         d*f2(j)*v0(j) - 2.d0*f3(j)*v0(j) + d*f3(j)*v0(j) - 
     $         2.d0*v1(j) + d*v1(j) + 2.d0*f1(j)*v1(j) - 
     $         d*f1(j)*v1(j) - 2.d0*f2(j)*v1(j) + d*f2(j)*v1(j) - 
     $         2.d0*f3(j)*v1(j) + 
     $         ((-2.d0 + d)*(-1.d0 + f1(j)) - 4.d0*f2(j))*v2(j))/2.d0 + 
     $      ((-3.d0 + d)*(ia23(j)*v0(j) + 
     $           ia22(j)*((-1.d0 + f1(j))*v0(j) + v1(j)) + 
     $           ia21(j)*((1.d0 - f1(j) + f2(j))*v0(j) + 
     $              (-1.d0 + f1(j))*v1(j) + v2(j))))/ia20(j)**2
      end do
      call inhom(v3, coeff1, coeff2, Delta, ny)


C      open(unit=10,file='testU_tayl.dat',status='unknown')
C            do j=1,ny
C                  write(10,*) u0(j), u1(j), u2(j), u3(j)
C            end do
C      close(10)
C      open(unit=10,file='testV_tayl.dat',status='unknown')
C            do j=1,ny
C                  write(10,*) v0(j), v1(j), v2(j), v3(j)
C            end do
C      close(10)
C      open(unit=10,file='testf_tayl.dat',status='unknown')
C            do j=1,ny
C                  write(10,*) f0(j), f1(j), f2(j), f3(j)
C            end do
C      close(10)
C      open(unit=10,file='testia2_tayl.dat',status='unknown')
C            do j=1,ny
C                  write(10,*) ia20(j), ia21(j), ia22(j), ia23(j)
C            end do
C      close(10)

C     Put together the expansion around x=xp=1.
      do j=1,ny
         f(j)  = f0(j) + (x-xp) * f1(j) + (x-xp)**2 * f2(j)
     $       + (x-xp)**3 * f3(j)
         ia2(j) = ia20(j) + (x-xp) * ia21(j) + (x-xp)**2 * ia22(j)
     $       + (x-xp)**3 * ia23(j)
         u(j)  = u0(j)+ (x-xp) * u1(j) + (x-xp)**2 * u2(j)
     $       + (x-xp)**3 * u3(j)
         v(j)  = v0(j)+ (x-xp) * v1(j) + (x-xp)**2 * v2(j)
     $       + (x-xp)**3 * v3(j)
      end do

      if (printtayl) then
            write(6,*) '***************************************'
            write(6,89) ' INFO: Taylor expansion at xright = ', x
            write(6,*) '***************************************'
 89         format (A,ES15.2)
            f10max = 0.d0
            f21max = 0.d0
            u10max = 0.d0
            u21max = 0.d0
            v10max = 0.d0
            v21max = 0.d0

            f0norm = 0.d0
            f1norm = 0.d0
            f2norm = 0.d0
            u0norm = 0.d0
            u1norm = 0.d0
            u2norm = 0.d0
            v0norm = 0.d0
            v1norm = 0.d0
            v2norm = 0.d0

            do j=1,ny
                  f10max = max(f10max,abs((x-xp)*f1(j)/f0(j)))
                  f21max = max(f21max,abs((x-xp)*f2(j)/f1(j)))
                  u10max = max(u10max,abs((x-xp)*u1(j)/u0(j)))
                  u21max = max(u21max,abs((x-xp)*u2(j)/u1(j)))
                  v10max = max(v10max,abs((x-xp)*v1(j)/v0(j)))
                  v21max = max(v21max,abs((x-xp)*v2(j)/v1(j)))

                  f0norm = f0norm + f0(j)**2
                  f1norm = f1norm + f1(j)**2
                  f2norm = f2norm + f2(j)**2
                  
                  u0norm = u0norm + u0(j)**2
                  u1norm = u1norm + u1(j)**2
                  u2norm = u2norm + u2(j)**2

                  v0norm = v0norm + v0(j)**2
                  v1norm = v1norm + v1(j)**2
                  v2norm = v2norm + v2(j)**2
                  
            end do

            f0norm = sqrt(f0norm)
            f1norm = sqrt(f1norm)
            f2norm = sqrt(f2norm)

            u0norm = sqrt(u0norm)
            u1norm = sqrt(u1norm)
            u2norm = sqrt(u2norm)
            v0norm = sqrt(v0norm)
            v1norm = sqrt(v1norm)
            v2norm = sqrt(v2norm)

            write(6,*) 'max((x-1) * f1 / f0) = ', f10max
            write(6,*) 'max((x-1) * f2 / f1) = ', f21max
            write(6,*)
            write(6,*) 'max((x-1) * u1 / u0) = ', u10max
            write(6,*) 'max((x-1) * u2 / u1) = ', u21max
            write(6,*)
            write(6,*) 'max((x-1) * v1 / v0) = ', v10max
            write(6,*) 'max((x-1) * v2 / v1) = ', v21max
            write(6,*)
            write(6,*) '((x-1) * f1norm) / f0norm = ', 
     $                  ((x-xp)*f1norm)/f0norm
            write(6,*) '((x-1) * f2norm) / f1norm = ', 
     $                  ((x-xp)*f2norm)/f1norm
            write(6,*)
            write(6,*) '((x-1) * u1norm) / u0norm = ', 
     $                  ((x-xp)*u1norm)/u0norm
            write(6,*) '((x-1) * u2norm) / u1norm = ', 
     $                  ((x-xp)*u2norm)/u1norm
            write(6,*)
            write(6,*) '((x-1) * v1norm) / v0norm = ', 
     $                  ((x-xp)*v1norm)/v0norm
            write(6,*) '((x-1) * v2norm) / v1norm = ', 
     $                  ((x-xp)*v2norm)/v1norm

      end if


      call yfromfields(ny, u, v, f, y)

      y(5) = Delta

C     Write out taylor coefficients      
      if (debug) then
C        possible debugging section            
         continue
      end if

      return
      end

