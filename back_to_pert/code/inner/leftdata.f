      subroutine leftdata(ny, d, x, Delta, fc, psic, y, debug, tstep)

      implicit none
      integer ny
      double precision x, Delta, d, fc(ny), psic(ny), y(ny)
      logical debug

      integer nymax, j, tstep
      include '../nymax.inc'
      double precision coeff1(nymax), coeff2(nymax),
     $     f(nymax), mu(nymax), pi(nymax), psi(nymax),
     $     f0(nymax),             pi0(nymax), psi0(nymax),
     $     f2(nymax), mu2(nymax), pi2(nymax), psi2(nymax),
     $     f4(nymax), mu4(nymax), pi4(nymax), psi4(nymax),
     $     df0dtau(nymax), dpi0dtau(nymax), dpsi0dtau(nymax), 
     $     d2f0dtau2(nymax), d2pi0dtau2(nymax), d2psi0dtau2(nymax), 
     $     d3f0dtau3(nymax), d3pi0dtau3(nymax), d3psi0dtau3(nymax),
     $     d4psi0dtau4(nymax), junkxi(nymax), dpi2dtau(nymax),
     $     u1(nymax), du1dtau(nymax), v2(nymax), u2(nymax), 
     $     v3(nymax), u3(nymax), v1(nymax), ia22(nymax), 
     $     u(nymax), v(nymax), d2u1dtau2(nymax), ia2(nymax)


      if (debug) write(6,*) 'Writing Taylor coeffs at x=', x


C *** Order 0 ***

      do j=1,ny
         f0(j) = fc(j)
         psi0(j) = psic(j)
      end do
     
    
C     Free data.
      call fourdiff1(ny, f0, df0dtau, Delta)
      call fourdiff1(ny, df0dtau, d2f0dtau2, Delta)
      call fourdiff1(ny, d2f0dtau2, d3f0dtau3, Delta)
      call fourdiff1(ny, psi0, dpsi0dtau, Delta)
      call fourdiff1(ny, dpsi0dtau, d2psi0dtau2, Delta)
      call fourdiff1(ny, d2psi0dtau2, d3psi0dtau3, Delta)
      call fourdiff1(ny, d3psi0dtau3, d4psi0dtau4, Delta)
      
     
      
C     Solve a linear inhomogeneous ODE for u1.
      do j=1,ny
         coeff1(j) = 1.d0
         coeff2(j) = -(d - 1.d0) * psi0(j) * f0(j)
      end do
      call inhom(u1, coeff1, coeff2, Delta, ny)
      call fourdiff1(ny, u1, du1dtau, Delta)
      call fourdiff1(ny, du1dtau, d2u1dtau2, Delta)

      do j=1,ny
         v1(j)=u1(j)
      end do
      
C *** Order 2 ***

      do j=1,ny
         f2(j) = ((d - 3.d0) * (d - 2.d0)**3 * f0(j) * u1(j)**2) 
     $           / (8.d0 * (d - 1.d0))
         u2(j) = - psi0(j)
         v2(j) = psi0(j)
         ia22(j) = -(d - 2.d0)**3*u1(j)**2/(4.d0*(d - 1.d0))
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
      


C     Put together the expansion around x=0. 
      do j=1,ny
         f(j) = f0(j)  + x**2 * f2(j)   
         u(j) = x * u1(j) + x**2 * u2(j)   
     $                    + x**3 * u3(j)
         v(j) = x * v1(j) + x**2 * v2(j)  
     $                    + x**3 * v3(j)
         ia2(j) = 1.d0 + x**2 * ia22(j)
      end do
      
      

      do j=1,ny
         junkxi(j) = 0.0d0
      end do

      call yfromfields(ny, u, v, f, y)
      y(5) = Delta
      
      
C     Output taylor coefficients
      if (debug) then
         if (tstep.eq.1) then
            open(unit=11,file='L_taylor/fleft.junk',status='new')
            open(unit=12,file='L_taylor/uleft.junk',status='new')
            open(unit=13,file='L_taylor/vleft.junk',status='new')
            open(unit=14,file='L_taylor/ia2left.junk',status='new')
      
               do j=1,ny
                  write(11,*) f0(j), f2(j)
                  write(12,*) u1(j), u2(j), u3(j)
                  write(13,*) v1(j), v2(j), v3(j)
                  write(14,*) ia22(j)
               end do
      
            do j=11,14
               close(j)
            end do
      
         else
         
            call changeresolution(ny, ny/2, f0, f0)
            call changeresolution(ny, ny/2, f2, f2)
            
            call changeresolution(ny, ny/2, u1, u1)
            call changeresolution(ny, ny/2, u2, u2)
            call changeresolution(ny, ny/2, u3, u3)
            
            call changeresolution(ny, ny/2, v1, v1)
            call changeresolution(ny, ny/2, v2, v2)
            call changeresolution(ny, ny/2, v3, v3)
           
            call changeresolution(ny, ny/2, ia22, ia22)
         
            open(unit=11,file='L_taylor/fleft.junk',status='new')
            open(unit=12,file='L_taylor/uleft.junk',status='new')
            open(unit=13,file='L_taylor/vleft.junk',status='new')
            open(unit=14,file='L_taylor/ia2left.junk',status='new')
      
               do j=1,ny/2
                  write(11,*) f0(j), f2(j)
                  write(12,*) u1(j), u2(j), u3(j)
                  write(13,*) v1(j), v2(j), v3(j)
                  write(14,*) ia22(j)
               end do
      
            do j=11,14
               close(j)
            end do
         end if
      end if
      


      return
      end
