      subroutine leftdatalin(ny, d, x, Delta, fc, psic, 
     $                      delfc, delpsic, y, lamb, debug)

      implicit none
      integer ny
      double precision x, Delta, lamb, d, fc(ny), psic(ny), y(ny),
     $     delfc(ny), delpsic(ny)
      logical debug

      integer nymax, j
      include '../nymax.inc'
      double precision coeff1(nymax), coeff2(nymax),
     $     f(nymax), u(nymax), v(nymax), ia2(nymax),
     $     f0(nymax), f2(nymax), 
     $     u1(nymax), u2(nymax), u3(nymax),
     $     v1(nymax), v2(nymax), v3(nymax),
     $     ia22(nymax), 
     $     psi0(nymax),
     $     df0dtau(nymax), d2f0dtau2(nymax),
     $     dpsi0dtau(nymax), 
     $     du1dtau(nymax), d2u1dtau2(nymax) 
           
C     Background variables
       double precision coeff1B(nymax), coeff2B(nymax),
     $     f0B(nymax), psi0B(nymax), u1B(nymax),
     $     du1dtauB(nymax), d2u1dtau2B(nymax),
     $     df0dtauB(nymax), d2f0dtau2B(nymax),
     $     dpsi0dtauB(nymax)


      if (debug) write(6,*) 'leftdata debug at x=', x

      

C *** Order 0 and background ***

      do j=1,ny
         f0(j) = delfc(j)
         psi0(j) = delpsic(j)
         f0B(j) = fc(j)
         psi0B(j) = psic(j)
      end do
     
    
C     Free data.
      call fourdiff1(ny, f0, df0dtau, Delta)
      call fourdiff1(ny, df0dtau, d2f0dtau2, Delta)
      call fourdiff1(ny, psi0, dpsi0dtau, Delta)
     
      call fourdiff1(ny, f0B, df0dtauB, Delta)
      call fourdiff1(ny, df0dtauB, d2f0dtau2B, Delta)
      call fourdiff1(ny, psi0B, dpsi0dtauB, Delta)
      
C     Solve a linear inhomogeneous ODE for u1 and u1B.
      do j=1,ny
         coeff1B(j) = 1.d0
         coeff2B(j) = -(d - 1.d0) * psi0B(j) * f0B(j)
         coeff1(j) = 1.d0 + lamb
         coeff2(j) = -(d - 1.d0) * (psi0(j) * f0B(j) + psi0B(j) * f0(j))
      end do
      call inhom(u1B, coeff1B, coeff2B, Delta, ny)
      call inhom(u1, coeff1, coeff2, Delta, ny)
      
      call fourdiff1(ny, u1, du1dtau, Delta)
      call fourdiff1(ny, du1dtau, d2u1dtau2, Delta)
      call fourdiff1(ny, u1B, du1dtauB, Delta)
      call fourdiff1(ny, du1dtauB, d2u1dtau2B, Delta)

      do j=1,ny
         v1(j)=u1(j)
      end do
C *** Order 2 ***

      do j=1,ny
         f2(j) = ((d - 3.d0) * (d - 2.d0)**3 * u1B(j) 
     $            * (2.d0*u1(j)*f0B(j) + f0(j)*u1B(j))) 
     $           / (8.d0 * (d - 1.d0))
         u2(j) = - psi0(j)
         v2(j) = psi0(j)
         ia22(j) = -((d - 2.d0)**3*u1(j)*u1B(j))/(2.d0*(d - 1.d0))
      end do


C *** Order 3 ***

      do j=1,ny
         u3(j)=(u1(j)*f0B(j)*(- 4.d0*(1.d0 + lamb)**2 - 3.d0*(d - 3.d0)
     $                          *(d - 2.d0)**3 * f0B(j)**2 * u1B(j)**2)
     $            + 4.d0*(d- 1.d0)*psi0B(j)*(f0(j)*((lamb - 3.d0)*f0B(j)
     $                                               + 2.d0*df0dtauB(j))
     $                                        -f0B(j)*df0dtau(j))
     $            + 4.d0*f0B(j)*((d - 1.d0)*psi0(j)*((3.d0 + 2.d0*lamb)
     $                           *f0B(j) - df0dtauB(j)) + d2u1dtau2(j))
     $            + 8.d0*f0(j)*(u1B(j) - d2u1dtau2B(j)))
     $           /( 8.d0 * (d - 1.d0) * (f0B(j))**3 ) 
         v3(j)=(u1(j)*f0B(j)*(- 4.d0*(1.d0 + lamb)**2 - 3.d0*(d - 3.d0)
     $                          *(d - 2.d0)**3 * f0B(j)**2 * u1B(j)**2)
     $            + 4.d0*(d- 1.d0)*psi0B(j)*(f0(j)*((lamb - 3.d0)*f0B(j)
     $                                               + 2.d0*df0dtauB(j))
     $                                        -f0B(j)*df0dtau(j))
     $            + 4.d0*f0B(j)*((d - 1.d0)*psi0(j)*((3.d0 + 2.d0*lamb)
     $                           *f0B(j) - df0dtauB(j)) + d2u1dtau2(j))
     $            + 8.d0*f0(j)*(u1B(j) - d2u1dtau2B(j)))
     $           /(8.d0 * (d - 1.d0) * f0B(j)**3) 
      end do
      

C     Put together the expansion around x=0. 
      do j=1,ny
         f(j) = f0(j)  + x**2 * f2(j)   
         u(j) = x * u1(j) + x**2 * u2(j)   
     $                    + x**3 * u3(j)
         v(j) = x * v1(j) + x**2 * v2(j)  
     $                    + x**3 * v3(j)
         ia2(j) = x**2 *ia22(j)
      end do
      
      
      
    
C     Debugging output.
      if (debug) then
      
      open(unit=10,file='leftpert.junk',status='new')
         do j=1,ny
            write(10,*) f(j), u(j), v(j), ia2(j)
         end do
      close(10)
      open(unit=11,file='fleftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) f0(j), f2(j)
         end do
      close(11)
      open(unit=11,file='uleftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) u1(j), u2(j), u3(j)
         end do
      close(11)
      open(unit=11,file='vleftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) v1(j), v2(j), v3(j)
         end do
      close(11)
      open(unit=11,file='ia2leftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) ia22(j)
         end do
      close(11)
      
      end if

      call yfromfields(ny, u, v, f, y)
    

      return
      end
