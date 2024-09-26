      subroutine rightdatalin(ny, d, x, Delta, Up, delUp, 
     $                        y, lamb, debug)

      implicit none
      integer ny
      double precision x, d, Delta, lamb, 
     $       Up(ny), delUp(ny), y(ny)
      logical debug

      integer nymax, j
      include '../nymax.inc'
      double precision xp
      double precision coeff1(nymax), coeff2(nymax),
     $     u(nymax), v(nymax), f(nymax), ia2(nymax),
     $     u0(nymax), v0(nymax), f0(nymax), ia20(nymax),
     $     u1(nymax), v1(nymax), f1(nymax), ia21(nymax),
     $     u2(nymax), v2(nymax), f2(nymax), ia22(nymax),
     $     du0dtau(nymax), d2u0dtau2(nymax)
 
C     Background valiables
       double precision u0B(nymax), v0B(nymax), f0B(nymax),
     $     v1B(nymax), v2B(nymax), ia20B(nymax),
     $     du0dtauB(nymax), d2u0dtau2B(nymax)

     
      if (debug) write(6,*) 'rightdata debug at x=', x

      xp = 1.d0

C ******************************************
C **** Load background variables ***********
C ******************************************

      do j=1,ny
         f0B(j) = xp
         u0B(j) = Up(j)
      end do
      call fourdiff1(ny, u0B, du0dtauB, Delta)
      call fourdiff1(ny, du0dtauB, d2u0dtau2B, Delta)   
      
      open(unit=10,file='taylor_coeff_R.dat',status='old')
         do j=1,ny
            read(10,*) ia20B(j), v0B(j), v1B(j), v2B(j)
         end do
      close(10)

C *** Order 0 ***

      do j=1,ny
         f0(j) = 0.d0
         u0(j) = delUp(j)
      end do
      call fourdiff1(ny, u0, du0dtau, Delta)
      call fourdiff1(ny, du0dtau, d2u0dtau2, Delta)   

C     Solve a linear inhomogeneous ODE for ia20 
      do j=1,ny
         coeff1(j) = 3.d0 - d + lamb 
     $               -((- 2.d0 + d)**3 * u0B(j)**2)/4.d0
         coeff2(j) = -((d - 2.d0)**3 
     $                 * ia20B(j)*u0(j)*u0B(j))/2.d0
      end do
      call inhom(ia20, coeff1, coeff2, Delta, ny)

C     Solve a linear inhomogeneous ODE for v0 
      do j=1,ny
         coeff1(j) = (6.d0 - 2.d0*d + 
     $               (- 2.d0 +d+ 2.d0*lamb)*ia20B(j))
     $               /(2.d0*ia20B(j))
         coeff2(j) = ((d - 2.d0)*u0(j))/2.d0 
     $              +((d - 3.d0)*ia20(j)*v0B(j))/ia20B(j)**2
      end do
      call inhom(v0, coeff1, coeff2, Delta, ny)

C *** Order 1 ***

C     f1, u1, ia21 are given algebraically
      do j=1,ny
         f1(j) = -((d - 3.d0)*ia20(j))/ia20B(j)**2  
         ia21(j) = -((d - 2.d0)**3 * ia20B(j)*(u0(j)*u0B(j) 
     $              + v0(j)*v0B(j)))/4.d0 - (ia20(j)*(8.d0*(d - 3.d0) 
     $              + (d - 2.d0)**3 * (u0B(j)**2 + v0B(j)**2)))/8.d0 
         u1(j) = ((d - 2.d0 + (6.d0 - 2.d0*d)/ia20B(j))*u0(j) 
     $            - 2.d0*(du0dtau(j) + lamb*u0(j)) 
     $            + (2.d0*(d - 3.d0)*ia20(j)*u0B(j))/ia20B(j)**2 
     $               + (d - 2.d0)*v0(j))/4.d0  
      end do

C     Solve a linear inhomogeneous ODE for v1. 
      do j=1,ny
         coeff1(j) = (3.d0*d)/2.d0 + lamb 
     $              + (6.d0 - 2.d0*d)/ia20B(j) - 3.d0
         coeff2(j) =  (4.d0*du0dtau(j)*ia20B(j)**2 - 
     $       2.d0*d*du0dtau(j)*ia20B(j)**2 - 12.d0*ia20(j)*u0B(j) + 
     $       10.d0*d*ia20(j)*u0B(j) - 2.d0*d**2 * ia20(j)*u0B(j) + 
     $       12.d0*d*ia20(j)*v0B(j) - 4.d0*d**2*ia20(j)*v0B(j) + 
     $       24.d0*ia20(j)*u0B(j)**2 * v0B(j) - 
     $       44.d0*d*ia20(j)*u0B(j)**2 * v0B(j) + 
     $       30.d0* d**2 * ia20(j)*u0B(j)**2 * v0B(j) - 
     $       9.d0* d**3 * ia20(j)*u0B(j)**2 * v0B(j) + 
     $       d**4 * ia20(j)*u0B(j)**2 * v0B(j) + 
     $       24.d0*ia20(j)*v0B(j)**3 - 44.d0*d*ia20(j)*v0B(j)**3 + 
     $       30.d0*d**2 * ia20(j)*v0B(j)**3 - 
     $       9.d0 * d**3 * ia20(j)*v0B(j)**3  
     $       + d**4 * ia20(j)*v0B(j)**3 -(d - 2.d0)*ia20B(j)*u0(j)
     $       *((- 6.d0 + 3.d0*d + 2.d0*lamb)*ia20B(j) + 
     $          2.d0*(d - 3.d0)*((d - 2.d0)**2*u0B(j)*v0B(j) - 1.d0)) 
     $       - ia20B(j)*v0(j)*
     $        (3.d0*(d - 2.d0)**2 * ia20B(j) + 
     $          (d - 3.d0)*(- 4.d0*d + (d - 2.d0)**3 * u0B(j)**2 + 
     $             3.d0*(d - 2.d0)**3 * v0B(j)**2)) - 
     $       48.d0*ia20(j)*v1B(j) + 16.d0*d*ia20(j)*v1B(j))
     $       /(8.d0*ia20B(j)**2)    
      end do
      call inhom(v1, coeff1, coeff2, Delta, ny)

C *** Order 2 ***

C     f2, u2, ia22 are given algebraically
      do j=1,ny
         f2(j) = -((d - 3.d0)*(d - 2.d0)*
     $           (- 2.d0*(d - 2.d0)**2*ia20B(j)*
     $           (u0(j)*u0B(j) + v0(j)*v0B(j)) + 
     $           ia20(j)*(- 8.d0 + (d - 2.d0)**2*u0B(j)**2 + 
     $           (d - 2.d0)**2*v0B(j)**2)))/(16.d0*ia20B(j)**2)
         ia22(j) = 4.d0*ia20B(j)*u0(j)*u0B(j) + 
     $        ((d - 2.d0)*ia20(j)*(64.d0*(d - 3.d0) + 
     $        (d - 2.d0)**2*((d - 2.d0)**3*u0B(j)**4 + 
     $           u0B(j)*(8.d0*du0dtauB(j) - 4.d0*(d - 2.d0)*v0B(j)) + 
     $           2.d0*u0B(j)**2 *
     $            (- 16.d0 + 6.d0*d + (d - 2.d0)**3 * v0B(j)**2) + 
     $           v0B(j)*(8.d0*(- 5.d0 + 2.d0*d)*v0B(j) + 
     $              (d - 2.d0)**3 * v0B(j)**3 - 16.d0*v1B(j)))) + 
     $     4.d0*(ia20B(j)*u0(j)*
     $         ((d - 2.d0)**6 * u0B(j)**3 + 
     $           (d - 2.d0)**3 * (2.d0*du0dtauB(j) 
     $            - (d - 2.d0)*v0B(j)) + 
     $        u0B(j)*(2.d0*d*(- 120.d0 + d*(84.d0 
     $                         + d*(3.d0*d - 26.d0))) + 
     $          + 2.d0*(d - 2.d0)**3 * lamb 
     $                   + (d - 2.d0)**6 * v0B(j)**2))
     $         + (d - 2.d0)**3 *
     $         (2.d0*ia20B(j)*(du0dtau(j)*u0B(j) - 2.d0*v0B(j)*v1(j)) 
     $         + v0(j)*(- 4.d0*(d - 3.d0)*v0B(j) + 
     $              ia20B(j)*
     $               (-((d - 2.d0)*u0B(j)) + 
     $                 (8.d0*d + (d - 2.d0)**3 * u0B(j)**2)*v0B(j) + 
     $                 (d - 2.d0)**3*v0B(j)**3 - 
     $                 4.d0*(5.d0*v0B(j) + v1B(j)))))))/128.d0
         u2(j) = (36.d0*du0dtauB(j)*ia20(j)*ia20B(j) - 
     $    12.d0*d*du0dtauB(j)*ia20(j)*ia20B(j) - 
     $    36.d0*du0dtau(j)*ia20B(j)**2 + 
     $    12.d0*d*du0dtau(j)*ia20B(j)**2 + 
     $    4.d0*d2u0dtau2(j)*ia20B(j)**3 + 
     $    24.d0*du0dtau(j)*ia20B(j)**3 - 
     $    8.d0*d*du0dtau(j)*ia20B(j)**3 + 
     $    8.d0*lamb*du0dtau(j)*ia20B(j)**3 - 216.d0*ia20(j)*u0B(j) + 
     $    144.d0*d*ia20(j)*u0B(j) - 24.d0*d**2*ia20(j)*u0B(j) + 
     $    120.d0*ia20(j)*ia20B(j)*u0B(j) - 
     $    70.d0*d*ia20(j)*ia20B(j)*u0B(j) + 
     $    10.d0*d**2*ia20(j)*ia20B(j)*u0B(j) + 
     $    48.d0*ia20(j)*ia20B(j)*u0B(j)**3 - 
     $    88.d0*d*ia20(j)*ia20B(j)*u0B(j)**3 + 
     $    60.d0*d**2*ia20(j)*ia20B(j)*u0B(j)**3 - 
     $    18.d0*d**3*ia20(j)*ia20B(j)*u0B(j)**3 + 
     $    2.d0*d**4*ia20(j)*ia20B(j)*u0B(j)**3 + 
     $    12.d0*ia20(j)*ia20B(j)*v0B(j) - 
     $    10.d0*d*ia20(j)*ia20B(j)*v0B(j) + 
     $    2.d0*d**2*ia20(j)*ia20B(j)*v0B(j) + 
     $    24.d0*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $    44.d0*d*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 + 
     $    30.d0*d**2*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $    9.d0*d**3*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 + 
     $    d**4*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $    2.d0*(d - 2.d0)*ia20B(j)**2*v0(j)*
     $     (d - 3.d0 + 2.d0*ia20B(j) + 
     $       (d - 3.d0)*(d - 2.d0)**2*u0B(j)*v0B(j)) - 
     $    ia20B(j)*u0(j)*(- 12.d0*(d - 3.d0)**2 + 
     $       4.d0*(d - 2.d0 + 2.d0*d*lamb - lamb*(6.d0 + lamb))*
     $        ia20B(j)**2 + 
     $       (d - 3.d0)*ia20B(j)*
     $        (2.d0*(- 20.d0 + 5.d0*d - 6.d0*lamb) + 
     $          6.d0*(d - 2.d0)**3*u0B(j)**2 + (d - 2.d0)**3*v0B(j)**2)
     $       ) - 8.d0*ia20B(j)**3*v1(j) + 4.d0*d*ia20B(j)**3*v1(j))
     $    /(32.d0*ia20B(j)**3)
      end do
      
      
C     Solve a linear inhomogeneous ODE for v2. 
      do j=1,ny
         coeff1(j) = - 5.d0 + (5.d0*d)/2.d0 
     $               + lamb + (9.d0 - 3.d0*d)/ia20B(j)
         
         coeff2(j) = (ia20B(j)*v0(j)*(16.d0*(d - 2.d0)**2*d*ia20B(j)**2 
     $       + 8.d0*(d - 3.d0)**2*(d - 2.d0)**3*
     $        (u0B(j)**2 + 4.d0*v0B(j)**2) - 
     $       (d - 3.d0)*ia20B(j)*
     $        (64.d0 + 16.d0*d**2 + (d - 2.d0)**6*u0B(j)**4 + 
     $          1152.d0*v0B(j)**2 - 1952.d0*d*v0B(j)**2 + 
     $          1200.d0*d**2*v0B(j)**2 - 312.d0*d**3*v0B(j)**2 + 
     $          28.d0*d**4*v0B(j)**2 + 576.d0*v0B(j)**4 - 
     $          1728.d0*d*v0B(j)**4 + 2160.d0*d**2*v0B(j)**4 - 
     $          1440.d0*d**3*v0B(j)**4 + 540.d0*d**4*v0B(j)**4 - 
     $          108.d0*d**5*v0B(j)**4 + 9.d0*d**6*v0B(j)**4 + 
     $          4.d0*(d - 2.d0)**3*u0B(j)*
     $           (- 2.d0*du0dtauB(j) + (d - 2.d0)*v0B(j)) + 
     $          2.d0*(d - 2.d0)**3*u0B(j)**2*
     $           (- 12.d0 + 5.d0*(d - 2.d0)**3*v0B(j)**2) - 
     $          384.d0*v0B(j)*v1B(j) + 576.d0*d*v0B(j)*v1B(j) - 
     $          288.d0*d**2*v0B(j)*v1B(j) + 48.d0*d**3*v0B(j)*v1B(j)))
     $      - 4.d0*(d - 2.d0)*ia20B(j)*u0(j)*
     $     (- 4.d0*(d - 2.d0)*(d + 2.d0*lamb)*ia20B(j)**2 - 
     $       8.d0*(d - 3.d0)**2*(- 1.d0 + (d - 2.d0)**2*u0B(j)*v0B(j)) + 
     $       (d - 3.d0)*ia20B(j)*
     $        (24.d0 - 4.d0*d + 8.d0*lamb - (d - 2.d0)**3*u0B(j)**2 + 
     $          2.d0*(d - 2.d0)**5*u0B(j)**3*v0B(j) + 8.d0*v0B(j)**2 - 
     $          12.d0*d*v0B(j)**2 + 6.d0*d**2*v0B(j)**2 - 
     $          d**3*v0B(j)**2 + 
     $          2.d0*(d - 2.d0)**2*u0B(j)*
     $           (4.d0*(d - 4.d0)*v0B(j) + (d - 2.d0)**3*v0B(j)**3 + 
     $             4.d0*v1B(j)))) - 
     $    8.d0*(- 24.d0*du0dtauB(j)*ia20(j)*ia20B(j) + 
     $       20.d0*d*du0dtauB(j)*ia20(j)*ia20B(j) - 
     $       4.d0*d**2*du0dtauB(j)*ia20(j)*ia20B(j) + 
     $       24.d0*du0dtau(j)*ia20B(j)**2 - 
     $       20.d0*d*du0dtau(j)*ia20B(j)**2 + 
     $       4.d0*d**2*du0dtau(j)*ia20B(j)**2 - 
     $       16.d0*du0dtau(j)*ia20B(j)**3 + 
     $       16.d0*d*du0dtau(j)*ia20B(j)**3 - 
     $       4.d0*d**2*du0dtau(j)*ia20B(j)**3 + 
     $       144.d0*ia20(j)*u0B(j) - 168.d0*d*ia20(j)*u0B(j) + 
     $       64.d0*d**2*ia20(j)*u0B(j) - 8.d0*d**3*ia20(j)*u0B(j) - 
     $       120.d0*ia20(j)*ia20B(j)*u0B(j) + 
     $       136.d0*d*ia20(j)*ia20B(j)*u0B(j) - 
     $       50.d0*d**2*ia20(j)*ia20B(j)*u0B(j) + 
     $       6.d0*d**3*ia20(j)*ia20B(j)*u0B(j) + 
     $       16.d0*f2(j)*ia20B(j)**3*u0B(j) - 
     $       8.d0*d*f2(j)*ia20B(j)**3*u0B(j) + 
     $       16.d0*ia20B(j)**3*u2(j) - 8.d0*d*ia20B(j)**3*u2(j) + 
     $       144.d0*ia20(j)*v0B(j) - 168.d0*d*ia20(j)*v0B(j) + 
     $       64.d0*d**2*ia20(j)*v0B(j) - 8.d0*d**3*ia20(j)*v0B(j) - 
     $       168.d0*ia20(j)*ia20B(j)*v0B(j) + 
     $       224.d0*d*ia20(j)*ia20B(j)*v0B(j) - 
     $       86.d0*d**2*ia20(j)*ia20B(j)*v0B(j) + 
     $       10.d0*d**3*ia20(j)*ia20B(j)*v0B(j) - 
     $       48.d0*f2(j)*ia20B(j)**2*v0B(j) + 
     $       16.d0*d*f2(j)*ia20B(j)**2*v0B(j) + 
     $       16.d0*f2(j)*ia20B(j)**3*v0B(j) - 
     $       8.d0*d*f2(j)*ia20B(j)**3*v0B(j) + 
     $       48.d0*ia20B(j)*ia22(j)*v0B(j) - 
     $       16.d0*d*ia20B(j)*ia22(j)*v0B(j) + 
     $       48.d0*du0dtauB(j)*ia20(j)*ia20B(j)*u0B(j)*v0B(j) - 
     $       88.d0*d*du0dtauB(j)*ia20(j)*ia20B(j)*u0B(j)*v0B(j) + 
     $       60.d0*d**2*du0dtauB(j)*ia20(j)*ia20B(j)*u0B(j)*
     $        v0B(j) - 18.d0*d**3*du0dtauB(j)*ia20(j)*ia20B(j)*
     $        u0B(j)*v0B(j) + 
     $       2.d0*d**4*du0dtauB(j)*ia20(j)*ia20B(j)*u0B(j)*
     $        v0B(j) - 216.d0*ia20(j)*u0B(j)**2*v0B(j) + 
     $       468.d0*d*ia20(j)*u0B(j)**2*v0B(j) - 
     $       402.d0*d**2*ia20(j)*u0B(j)**2*v0B(j) + 
     $       171.d0*d**3*ia20(j)*u0B(j)**2*v0B(j) - 
     $       36.d0*d**4*ia20(j)*u0B(j)**2*v0B(j) + 
     $       3.d0*d**5*ia20(j)*u0B(j)**2*v0B(j) + 
     $       24.d0*d*ia20(j)*ia20B(j)*u0B(j)**2*v0B(j) - 
     $       44.d0*d**2*ia20(j)*ia20B(j)*u0B(j)**2*v0B(j) + 
     $       30.d0*d**3*ia20(j)*ia20B(j)*u0B(j)**2*v0B(j) - 
     $       9.d0*d**4*ia20(j)*ia20B(j)*u0B(j)**2*v0B(j) + 
     $       d**5*ia20(j)*ia20B(j)*u0B(j)**2*v0B(j) + 
     $       48.d0*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $       112.d0*d*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 + 
     $       104.d0*d**2*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $       48.d0*d**3*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 + 
     $       11.d0*d**4*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $       d**5*ia20(j)*ia20B(j)*u0B(j)*v0B(j)**2 - 
     $       72.d0*ia20(j)*v0B(j)**3 + 156.d0*d*ia20(j)*v0B(j)**3 - 
     $       134.d0*d**2*ia20(j)*v0B(j)**3 + 
     $       57.d0*d**3*ia20(j)*v0B(j)**3 - 
     $       12.d0*d**4*ia20(j)*v0B(j)**3 + 
     $       d**5*ia20(j)*v0B(j)**3 - 
     $       48.d0*ia20(j)*ia20B(j)*v0B(j)**3 + 
     $       136.d0*d*ia20(j)*ia20B(j)*v0B(j)**3 - 
     $       148.d0*d**2*ia20(j)*ia20B(j)*v0B(j)**3 + 
     $       78.d0*d**3*ia20(j)*ia20B(j)*v0B(j)**3 - 
     $       20.d0*d**4*ia20(j)*ia20B(j)*v0B(j)**3 + 
     $       2.d0*d**5*ia20(j)*ia20B(j)*v0B(j)**3 + 
     $       ia20B(j)**2*(8.d0*(10.d0 - 9.d0*d + 2.d0*d**2)*ia20B(j) + 
     $          (d - 3.d0)*(- 16.d0*(d - 1.d0) + 
     $             3.d0*(d - 2.d0)**3*u0B(j)**2 + 
     $             3.d0*(d - 2.d0)**3*v0B(j)**2))*v1(j) - 
     $       24.d0*d*ia20(j)*ia20B(j)*v1B(j) + 
     $       8.d0*d**2*ia20(j)*ia20B(j)*v1B(j) + 
     $       16.d0*f2(j)*ia20B(j)**3*v1B(j) - 
     $       48.d0*ia20(j)*ia20B(j)*u0B(j)**2*v1B(j) + 
     $       88.d0*d*ia20(j)*ia20B(j)*u0B(j)**2*v1B(j) - 
     $       60.d0*d**2*ia20(j)*ia20B(j)*u0B(j)**2*v1B(j) + 
     $       18.d0*d**3*ia20(j)*ia20B(j)*u0B(j)**2*v1B(j) - 
     $       2.d0*d**4*ia20(j)*ia20B(j)*u0B(j)**2*v1B(j) - 
     $       144.d0*ia20(j)*ia20B(j)*v0B(j)**2*v1B(j) + 
     $       264.d0*d*ia20(j)*ia20B(j)*v0B(j)**2*v1B(j) - 
     $       180.d0*d**2*ia20(j)*ia20B(j)*v0B(j)**2*v1B(j) + 
     $       54.d0*d**3*ia20(j)*ia20B(j)*v0B(j)**2*v1B(j) - 
     $       6.d0*d**4*ia20(j)*ia20B(j)*v0B(j)**2*v1B(j) + 
     $       144.d0*ia20(j)*ia20B(j)*v2B(j) - 
     $       48.d0*d*ia20(j)*ia20B(j)*v2B(j)))/(128.d0*ia20B(j)**3)
      end do
      call inhom(v2, coeff1, coeff2, Delta, ny)

C     Put together the expansion around x=xp=1.
      do j=1,ny
         f(j)  = f0(j) + (x-xp) * f1(j) + (x-xp)**2 * f2(j)
         ia2(j) = ia20(j) + (x-xp) * ia21(j) + (x-xp)**2 * ia22(j)
         u(j)  = u0(j)+ (x-xp) * u1(j) + (x-xp)**2 * u2(j)
         v(j)  = v0(j)+ (x-xp) * v1(j) + (x-xp)**2 * v2(j)
      end do
      
      if (debug) then 
      open(unit=10,file='rightpert.junk',status='new')
         do j=1,ny
            write(10,*) f(j), u(j), v(j), ia2(j)
         end do
      close(10)
      open(unit=11,file='frightcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) f0(j), f1(j), f2(j)
         end do
      close(11)
      open(unit=11,file='urightcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) u0(j), u1(j), u2(j)
         end do
      close(11)
      open(unit=11,file='vrightcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) v0(j), v1(j), v2(j)
         end do
      close(11)
      open(unit=11,file='ia2rightcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) ia20(j), ia21(j), ia22(j)
         end do
      close(11)
      stop
      end if              

      call yfromfields(ny, u, v, f, y)

      return
      end

