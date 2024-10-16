      subroutine rightdata(ny, d, x, Delta, Up, y, debug, tstep)

      implicit none
      integer ny
      double precision x, d, Delta, Up(ny), y(ny)
      logical debug

      integer nymax, j, tstep
      include '../nymax.inc'
      double precision xp
      double precision coeff1(nymax), coeff2(nymax),
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     u0(nymax), v0(nymax), f0(nymax), ia20(nymax),
     $     u1(nymax), v1(nymax), f1(nymax), ia21(nymax),
     $     u2(nymax), v2(nymax), f2(nymax), ia22(nymax),
     $     du0dtau(nymax), d2u0dtau2(nymax)

      if (debug) write(6,*) 'Writing Taylor coeffs at x=', x

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
      end do
      
     
C     Write taylor coefficients to files, needed for linearized theory
      if (tstep.eq.1) then
         open(unit=10, file='bg_data/taylor_coeff_R.dat', status='new')
             do j=1,ny
                write(10,*) ia20(j), v0(j), v1(j), v2(j)
             end do
         close(10)
      else
         open(unit=10, file='bg_data/taylor_coeff_R.dat', status='new')
             call changeresolution(ny, ny/2, ia20, ia20)
             call changeresolution(ny, ny/2, v0, v0)
             call changeresolution(ny, ny/2, v1, v1)
             call changeresolution(ny, ny/2, v2, v2)
             do j=1,ny/2
                write(10,*) ia20(j), v0(j), v1(j), v2(j)
             end do
         close(10)
      end if

      call yfromfields(ny, u, v, f, y)

      y(5) = Delta
      
      
C     Write out taylor coefficients      
      if (debug) then
         call system('rm -rf R_taylor; mkdir R_taylor')
         if (tstep.eq.1) then
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
      
         else
         
            call changeresolution(ny, ny/2, f0, f0)
            call changeresolution(ny, ny/2, f1, f1)
            call changeresolution(ny, ny/2, f2, f2)
            
            call changeresolution(ny, ny/2, u0, u0)
            call changeresolution(ny, ny/2, u1, u1)
            call changeresolution(ny, ny/2, u2, u2)
           
            call changeresolution(ny, ny/2, ia21, ia21)
            call changeresolution(ny, ny/2, ia22, ia22)
         
            open(unit=11,file='R_taylor/fright.junk',status='new')
            open(unit=12,file='R_taylor/uright.junk',status='new')
            open(unit=13,file='R_taylor/vright.junk',status='new')
            open(unit=14,file='R_taylor/ia2right.junk',status='new')
      
               do j=1,ny/2
                  write(11,*) f0(j), f1(j), f2(j)
                  write(12,*) u0(j), u1(j), u2(j)
                  write(13,*) v0(j), v1(j), v2(j)
                  write(14,*) ia20(j), ia21(j), ia22(j)
               end do
      
            do j=11,14
               close(j)
            end do
         end if
      end if

      return
      end

