      subroutine leftdata(ny, d, x, Delta, fc, psic, y, debug)

      implicit none
      integer ny
      double precision x, Delta, d, fc(ny), psic(ny), y(ny)
      logical debug

      integer nymax, j
      include '../nymax.inc'
      double precision coeff1(nymax), coeff2(nymax),
     $     f(nymax), 
     $     f0(nymax), psi0(nymax),
     $     f2(nymax), 
     $     f4(nymax), 
     $     df0dtau(nymax), du1dtau(nymax), d2u1dtau2(nymax),
     $     d2f0dtau2(nymax), d3u1dtau3(nymax),
     $     d3f0dtau3(nymax), d4u1dtau4(nymax), 
     $     u1(nymax), v2(nymax), u2(nymax), 
     $     v3(nymax), u3(nymax), v1(nymax),
     $     u(nymax), v(nymax), u4(nymax),
     $     v4(nymax), u5(nymax), v5(nymax)

C     Taylor variables     
      double precision f20max, f42max,
     $     f0norm, f2norm, f4norm,
     $     u21max, u32max, u43max, u54max, 
     $     u1norm, u2norm, u3norm, u4norm, u5norm 
            
      if (debug) write(6,*) 'leftdata debug at x=', x


C *** Order 0 ***

      do j=1,ny
         f0(j) = fc(j)
         psi0(j) = psic(j)
      end do
     
    
C     Free data.
      call fourdiff1(ny, f0, df0dtau, Delta)
      call fourdiff1(ny, df0dtau, d2f0dtau2, Delta)
      call fourdiff1(ny, d2f0dtau2, d3f0dtau3, Delta)
      
      
C     Solve a linear inhomogeneous ODE for u1.
      do j=1,ny
         coeff1(j) = 1.d0
         coeff2(j) = -(d - 1.d0) * psi0(j) * f0(j)
      end do
      call inhom(u1, coeff1, coeff2, Delta, ny)
      call fourdiff1(ny, u1, du1dtau, Delta)
      call fourdiff1(ny, du1dtau, d2u1dtau2, Delta)
      call fourdiff1(ny, d2u1dtau2, d3u1dtau3, Delta)
      call fourdiff1(ny, d3u1dtau3, d4u1dtau4, Delta)

      do j=1,ny
         v1(j)=u1(j)
      end do
C *** Order 2 ***

      do j=1,ny
         f2(j) = ((d - 3.d0) * (d - 2.d0)**3 * f0(j) * u1(j)**2) 
     $           / (8.d0 * (d - 1.d0))
         u2(j) = - psi0(j)
         v2(j) = psi0(j)
      end do


C *** Order 3 ***

      do j=1,ny
         u3(j) = ((d - 3.d0) * (d - 2.d0)**3 * f0(j)**3 * u1(j)**3  
     $             + 4.d0*df0dtau(j) * (u1(j) + du1dtau(j))
     $             - 4.d0*f0(j) * (2.d0*u1(j) 
     $                             + 3.d0*du1dtau(j) + d2u1dtau2(j)))
     $             / (8.d0*(1.d0 - d) * f0(j)**3)
         v3(j) = ((d - 3.d0) * (d - 2.d0)**3 * f0(j)**3 * u1(j)**3  
     $             + 4.d0*df0dtau(j) * (u1(j) + du1dtau(j))
     $             - 4.d0*f0(j) * (2.d0*u1(j) 
     $                             + 3.d0*du1dtau(j) + d2u1dtau2(j)))
     $             / (8.d0*(1.d0 - d) * f0(j)**3)
      end do

C *** Order 4 ***

      do j=1,ny
            f4(j) = - ((d - 3.d0*d)*(d - 2.d0)**3*
     $          ((d - 2.d0)**3*(5.d0 - 6.d0*d + d**2)*f0(j)**3*u1(j)**4
     $              + 8*(d - 1.d0)*df0dtau(j)*u1(j)*
     $               (du1dtau(j) + u1(j)) - 
     $              8.d0*f0(j)*(du1dtau(j)**2 + 
     $                 ((d - 1.d0)*d2u1dtau2(j) + 
     $                    (-1.d0 + 3.d0*d)*du1dtau(j))*u1(j) + 
     $                 (-1.d0 + 2.d0*d)*u1(j)**2)))/
     $          (128.d0*(-1.d0 + d)**2*(1.d0 + d)*f0(j)**2)
            u4(j) = (-192.d0*(-1.d0 + d)*df0dtau(j)**2*du1dtau(j)
     $             *f0(j)**2 + 
     $           192.d0*(-1.d0 + d)*d2u1dtau2(j)*df0dtau(j)*f0(j)**3 + 
     $           64.d0*(-1.d0 + d)*d2f0dtau2(j)*du1dtau(j)*f0(j)**3 + 
     $           640.d0*(-1.d0 + d)*df0dtau(j)*du1dtau(j)*f0(j)**3 - 
     $           64.d0*(-1.d0 + d)*(6.d0*d2u1dtau2(j) + d3u1dtau3(j))*
     $            f0(j)**4 - 704.d0*(-1.d0 + d)*du1dtau(j)*f0(j)**4 - 
     $           192.d0*(-1.d0 + d)*df0dtau(j)**2*f0(j)**2*u1(j) + 
     $           64.d0*(-1.d0 + d)*d2f0dtau2(j)*f0(j)**3*u1(j) + 
     $           448.d0*(-1.d0 + d)*df0dtau(j)*f0(j)**3*u1(j) - 
     $           384.d0*(-1.d0 + d)*f0(j)**4*u1(j) + 
     $           32.d0*(-2.d0 + d)**3*(3.d0 - 7.d0*d + 2.d0*d**2)*
     $             f0(j)**6*u1(j)**2*
     $            (du1dtau(j) + u1(j)))/
     $         (128.d0*(-1.d0 + d)**2*(1.d0 + d)*f0(j)**7)

            v4(j) = - u4(j)
      end do

C *** Order 5 ***

      do j=1,ny
      u5(j) = (-64.d0*(-1.d0 + d)*d2f0dtau2(j)*d2u1dtau2(j)*f0(j)**2 - 
     $        16.d0*(-1.d0 + d)*(15.d0*d2f0dtau2(j) + d3f0dtau3(j))*
     $            du1dtau(j)*f0(j)**2 - 
     $           32.d0*(-1.d0 + d)*df0dtau(j)*
     $     (20.d0*d2u1dtau2(j) + 3.d0*d3u1dtau3(j) + 40.d0*du1dtau(j))*
     $            f0(j)**2 + 560.d0*(-1.d0 + d)*d2u1dtau2(j)*f0(j)**3 + 
     $           160.d0*(-1.d0 + d)*d3u1dtau3(j)*f0(j)**3 + 
     $           16.d0*(-1.d0 + d)*d4u1dtau4(j)*f0(j)**3 + 
     $           800.d0*(-1.d0 + d)*du1dtau(j)*f0(j)**3 - 
     $           16.d0*(-1.d0 + d)*(11.d0*d2f0dtau2(j) + d3f0dtau3(j))*
     $            f0(j)**2*u1(j) - 
     $           736.d0*(-1.d0 + d)*df0dtau(j)*f0(j)**2*u1(j) + 
     $           384.d0*(-1.d0 + d)*f0(j)**3*u1(j) - 
     $   8.d0*(-2.d0 + d)**3*(3.d0 - 13.d0*d + 4.d0*d**2)*du1dtau(j)**2*
     $            f0(j)**5*u1(j) + 
     $      8.d0*(-2.d0 + d)**3*(3.d0 - 13.d0*d + 4.d0*d**2)*df0dtau(j)*
     $            du1dtau(j)*f0(j)**4*u1(j)**2 - 
     $           8.d0*(-2.d0 + d)**3*(3.d0 - 13.d0*d + 4.d0*d**2)*
     $            (d2u1dtau2(j) + 5.d0*du1dtau(j))*f0(j)**5*u1(j)**2 + 
     $      8.d0*(-2.d0 + d)**3*(3.d0 - 13.d0*d + 4.d0*d**2)*df0dtau(j)*
     $            f0(j)**4*u1(j)**3 - 
     $       24.d0*(-2.d0 + d)**3*(3.d0 - 13.d0*d + 4.d0*d**2)*f0(j)**5*
     $            u1(j)**3 + (-2.d0 + d)**6*
     $     (3.d0 + 29.d0*d - 19.d0*d**2 + 3.d0*d**3)*f0(j)**7*u1(j)**5 - 
     $           240.d0*(-1.d0 + d)*df0dtau(j)**3*(du1dtau(j) + u1(j)) + 
     $           80.d0*(-1.d0 + d)*df0dtau(j)*f0(j)*
     $            (2.d0*d2f0dtau2(j)*du1dtau(j) + 
     $              df0dtau(j)*(3.d0*d2u1dtau2(j) + 11.d0*du1dtau(j)) + 
     $              2.d0*(d2f0dtau2(j) + 4.d0*df0dtau(j))*u1(j)))/
     $         (128.d0*(-1.d0 + d)**2*(1.d0 + d)*f0(j)**7)

      v5(j) = u5(j)
      end do


C     Put together the expansion around x=0. 
      do j=1,ny
         f(j) = f0(j)  + x**2 * f2(j) + x**4 * f4(j)  
         u(j) = x * u1(j) + x**2 * u2(j)   
     $          + x**3 * u3(j)
     $          + x**4 * u4(j) + x**5 * u5(j)
         v(j) = x * v1(j) + x**2 * v2(j)  
     $          + x**3 * v3(j)
     $          + x**4 * v4(j) + x**5 * v5(j)
      end do

      write(6,89) ' INFO: Taylor expansion at xleft = ', x
 89   format (A,ES15.2)
      f20max = 0.d0
      f42max = 0.d0
      u21max = 0.d0
      u32max = 0.d0
      u43max = 0.d0
      u54max = 0.d0

      f0norm = 0.d0
      f2norm = 0.d0
      f4norm = 0.d0
      u1norm = 0.d0
      u2norm = 0.d0
      u3norm = 0.d0
      u4norm = 0.d0
      u5norm = 0.d0

      do j=1,ny
            f20max = max(f20max,abs(x**2*f2(j)/f0(j)))
            f42max = max(f42max,abs(x**2*f4(j)/f2(j)))
            u21max = max(u21max,abs(x*u2(j)/u1(j)))
            u32max = max(u32max,abs(x*u3(j)/u2(j)))
            u43max = max(u43max,abs(x*u4(j)/u3(j)))
            u54max = max(u54max,abs(x*u5(j)/u4(j)))

            f0norm = f0norm + f0(j)**2
            f2norm = f2norm + f2(j)**2
            f4norm = f4norm + f4(j)**2
            
            u1norm = u1norm + u1(j)**2
            u2norm = u2norm + u2(j)**2
            u3norm = u3norm + u3(j)**2
            u4norm = u4norm + u4(j)**2
            u5norm = u5norm + u5(j)**2
      end do

      f0norm = sqrt(f0norm)
      f2norm = sqrt(f2norm)
      f4norm = sqrt(f4norm)

      u1norm = sqrt(u1norm)
      u2norm = sqrt(u2norm)
      u3norm = sqrt(u3norm)
      u4norm = sqrt(u4norm)
      u5norm = sqrt(u5norm)

      write(6,*) 'max(x^2 f2 / f0) = ', f20max
      write(6,*) 'max(x^2 f4 / f2) = ', f42max
      write(6,*)
      write(6,*) 'max(x * u2 / u1) = ', u21max
      write(6,*) 'max(x * u3 / u2) = ', u32max
      write(6,*) 'max(x * u4 / u3) = ', u43max
      write(6,*) 'max(x * u5 / u4) = ', u54max
      write(6,*)
      write(6,*) '(x^2 f2norm) / f0norm = ', (x**2*f2norm)/f0norm
      write(6,*) '(x^2 f4norm) / f2norm = ', (x**2*f4norm)/f2norm
      write(6,*)
      write(6,*) '(x * u2norm) / u1norm = ', (x*u2norm)/u1norm
      write(6,*) '(x * u3norm) / u2norm = ', (x*u3norm)/u2norm
      write(6,*) '(x * u4norm) / u3norm = ', (x*u4norm)/u3norm
      write(6,*) '(x * u5norm) / u4norm = ', (x*u5norm)/u4norm

      call yfromfields(ny, u, v, f, y)
      y(5) = Delta

C     Output taylor coefficients
      if (debug) then
        
            open(unit=11,file='L_taylor/fleft.junk',status='new')
            open(unit=12,file='L_taylor/uleft.junk',status='new')
            open(unit=13,file='L_taylor/vleft.junk',status='new')
      
               do j=1,ny
                  write(11,*) f0(j), f2(j)
                  write(12,*) u1(j), u2(j), u3(j)
                  write(13,*) v1(j), v2(j), v3(j)
               end do
      
            do j=11,13
               close(j)
            end do
      end if

      return
      end
