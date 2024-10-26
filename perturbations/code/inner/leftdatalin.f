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
     $     f0(nymax), f2(nymax), f4(nymax),
     $     u1(nymax), u2(nymax), u3(nymax), u4(nymax),
     $     v1(nymax), v2(nymax), v3(nymax), v4(nymax),
     $     u5(nymax), v5(nymax), ia22(nymax), ia24(nymax),
     $     psi0(nymax),
     $     df0dtau(nymax), d2f0dtau2(nymax),
     $     d3f0dtau3(nymax), 
     $     du1dtau(nymax), d2u1dtau2(nymax),
     $     d3u1dtau3(nymax), d4u1dtau4(nymax)
           
C     Background variables
       double precision coeff1B(nymax), coeff2B(nymax),
     $     f0B(nymax), psi0B(nymax), u1B(nymax),
     $     du1dtauB(nymax), d2u1dtau2B(nymax),
     $     d3u1dtau3B(nymax), d4u1dtau4B(nymax),
     $     df0dtauB(nymax), d2f0dtau2B(nymax),
     $     d3f0dtau3B(nymax)


C      if (debug) write(6,*) 'leftdata debug at x=', x

      

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
      call fourdiff1(ny, d2f0dtau2, d3f0dtau3, Delta)
     
      call fourdiff1(ny, f0B, df0dtauB, Delta)
      call fourdiff1(ny, df0dtauB, d2f0dtau2B, Delta)
      call fourdiff1(ny, d2f0dtau2B, d3f0dtau3B, Delta)
      
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
      call fourdiff1(ny, d2u1dtau2, d3u1dtau3, Delta)
      call fourdiff1(ny, d3u1dtau3, d4u1dtau4, Delta)
      call fourdiff1(ny, u1B, du1dtauB, Delta)
      call fourdiff1(ny, du1dtauB, d2u1dtau2B, Delta)
      call fourdiff1(ny, d2u1dtau2B, d3u1dtau3B, Delta)
      call fourdiff1(ny, d3u1dtau3B, d4u1dtau4B, Delta)

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
      
      
C *** Order 4 ***

      do j=1,ny
         f4(j) = -((-3.d0 + d)*(-2.d0 + d)**3*
     $         (8.d0*f0B(j)*(-2.d0*du1dtau(j)*du1dtauB(j)*f0B(j) + 
     $              ((-1.d0 + d)*
     $                  (df0dtauB(j)*du1dtau(j) + 
     $                    df0dtau(j)*du1dtauB(j)) + 
     $                 (-((-1.d0 + d)*d2u1dtau2(j)) + 
     $                    (1.d0 - 3.d0*d + 2.d0*lamb - 2.d0*d*lamb)*
     $                     du1dtau(j))*f0B(j))*u1B(j) + 
     $              (-1.d0 + d)*df0dtau(j)*u1B(j)**2) + 
     $           4.d0*f0B(j)*u1(j)*
     $            ((-5.d0 + d)*(-2.d0 + d)**3*(-1.d0 + d)*f0B(j)**3*
     $               u1B(j)**3 + 
     $              2.d0*(-1.d0 + d)*df0dtauB(j)*
     $               (du1dtauB(j) + (2.d0 + lamb)*u1B(j)) + 
     $              2.d0*f0B(j)*
     $               (-((-1.d0 + d)*d2u1dtau2B(j)) + 
     $                 (1.d0 - 3.d0*d - 2.d0*lamb)*du1dtauB(j) + 
     $                 (2.d0 + lamb + lamb**2 - 
     $                    d*(4.d0 + lamb*(3.d0 + lamb)))*u1B(j))) + 
     $         f0(j)*((-5.d0 + d)*(-2.d0 + d)**3*(-1.d0 + d)*f0B(j)**3*
     $               u1B(j)**4 - 
     $              16.d0*(-1.d0 + d)*df0dtauB(j)*u1B(j)*
     $               (du1dtauB(j) + u1B(j)) + 
     $              8.d0*f0B(j)*
     $               (du1dtauB(j)**2 + 
     $                 ((-1.d0 + d)*d2u1dtau2B(j) + 
     $                    (-1.d0 + 3*d + (-1.d0 + d)*lamb)*du1dtauB(j)
     $                    )*u1B(j) + 
     $             (-1.d0 + 2.d0*d + (-1.d0 + d)*lamb)*u1B(j)**2))))/
     $       (128.d0*(-1.d0 + d)**2*(1.d0 + d)*f0B(j)**3)


         u4(j) = (f0B(j)*u1(j)*(-6.d0*(-1.d0 + d)*(1.d0 + lamb)*
     $            df0dtauB(j)**2 + 
     $           2.d0*(-1.d0 + d)*(1.d0 + lamb)*
     $            (d2f0dtau2B(j) + (7.d0 + 3.d0*lamb)*df0dtauB(j))*
     $            f0B(j) - 2.d0*(-1.d0 + d)*
     $          (6.d0 + 11.d0*lamb + 6.d0*lamb**2 + lamb**3)*f0B(j)**2
     $          + (-2.d0 + d)**3*(3.d0 - 7.d0*d + 2.d0*d**2)*f0B(j)**4*
     $            u1B(j)*(2.d0*du1dtauB(j) + (3.d0 + lamb)*u1B(j))) + 
     $        f0B(j)*(-2.d0*(-1.d0 + d)*
     $            (3.d0*(2.d0 + lamb)*d2u1dtau2(j) + d3u1dtau3(j) + 
     $              (11.d0 + 12.d0*lamb + 3.d0*lamb**2)*du1dtau(j))*
     $            f0B(j)**2 + 
     $           (-2.d0 + d)**3*(3.d0 - 7.d0*d + 2.d0*d**2)*du1dtau(j)*
     $            f0B(j)**4*u1B(j)**2 - 
     $           6.d0*(-1.d0 + d)*df0dtauB(j)*
     $            (df0dtauB(j)*du1dtau(j) + 
     $              2.d0*df0dtau(j)*du1dtauB(j) + 
     $              2.d0*df0dtau(j)*u1B(j)) + 
     $           2.d0*(-1.d0 + d)*f0B(j)*
     $            (3.d0*d2u1dtau2B(j)*df0dtau(j) + 
     $              3.d0*d2u1dtau2(j)*df0dtauB(j) + 
     $              (d2f0dtau2B(j) + 
     $                 2.d0*(5.d0 + 3.d0*lamb)*df0dtauB(j))*du1dtau(j) + 
     $              d2f0dtau2(j)*du1dtauB(j) + 
     $              10.d0*df0dtau(j)*du1dtauB(j) + 
     $              2.d0*lamb*df0dtau(j)*du1dtauB(j) + 
     $              (d2f0dtau2(j) + (7.d0 + 2.d0*lamb)*df0dtau(j))*
     $               u1B(j))) + 
     $        f0(j)*(30.d0*(-1.d0 + d)*df0dtauB(j)**2*
     $            (du1dtauB(j) + u1B(j)) - 
     $           (-2.d0 + d)**3*(3.d0 - 7.d0*d + 2.d0*d**2)*f0B(j)**4*
     $            u1B(j)**2*(du1dtauB(j) + u1B(j)) + 
     $           2.d0*(-1.d0 + d)*f0B(j)**2*
     $            (3.d0*((6.d0 + lamb)*d2u1dtau2B(j) + 
     $                 d3u1dtau3B(j)) + 
     $              (33.d0 + 10.d0*lamb + lamb**2)*du1dtauB(j) + 
     $              (18.d0 + 7.d0*lamb + lamb**2)*u1B(j)) - 
     $           4.d0*(-1.d0 + d)*f0B(j)*
     $            (2.d0*d2f0dtau2B(j)*du1dtauB(j) + 
     $              df0dtauB(j)*
     $               (6.d0*d2u1dtau2B(j) + 
     $                 (20.d0 + 3.d0*lamb)*du1dtauB(j)) + 
     $              (2.d0*d2f0dtau2B(j) + 
     $                 (14.d0 + 3.d0*lamb)*df0dtauB(j))*u1B(j))))/
     $      (4.d0*(-1.d0 + d)**2*(1.d0 + d)*f0B(j)**6)

         v4(j) = - u4(j)

         ia24(j) = ((-2.d0 + d)**3*(f0B(j)*
     $           (-2.d0*du1dtau(j)*du1dtauB(j)*f0B(j) + 
     $             ((-1.d0 + d)*
     $                 (df0dtauB(j)*du1dtau(j) + 
     $                   df0dtau(j)*du1dtauB(j)) + 
     $                (-((-1.d0 + d)*d2u1dtau2(j)) + 
     $                   (1.d0 - 3.d0*d + 2.d0*lamb - 2.d0*d*lamb)*
     $                    du1dtau(j))*f0B(j))*u1B(j) + 
     $             (-1.d0 + d)*df0dtau(j)*u1B(j)**2) + 
     $          f0B(j)*u1(j)*
     $           ((-2.d0 + d)**4*(-1.d0 + d)*f0B(j)**3*u1B(j)**3 + 
     $             (-1.d0 + d)*df0dtauB(j)*
     $              (du1dtauB(j) + (2.d0 + lamb)*u1B(j)) + 
     $             f0B(j)*(-((-1.d0 + d)*d2u1dtau2B(j)) + 
     $                (1.d0 - 3.d0*d - 2.d0*lamb)*du1dtauB(j) + 
     $                (2.d0 + lamb + lamb**2 - 
     $                   d*(4.d0 + lamb*(3.d0 + lamb)))*u1B(j))) + 
     $          f0(j)*(-3.d0*(-1.d0 + d)*df0dtauB(j)*u1B(j)*
     $              (du1dtauB(j) + u1B(j)) + 
     $             f0B(j)*(2.d0*du1dtauB(j)**2 + 
     $                (2.d0*(-1.d0 + d)*d2u1dtau2B(j) + 
     $                   (-2.d0 - lamb + d*(6.d0 + lamb))*du1dtauB(j))
     $                  *u1B(j) + 
     $                (-2.d0 - lamb + d*(4.d0 + lamb))*u1B(j)**2))))/
     $      (4.d0*(-1.d0 + d)**2*(1.d0 + d)*f0B(j)**4)

      end do

C *** Order 5 ***
      
      do j=1,ny
         u5(j) = (f0B(j)*u1(j)*(-240.d0*(-1.d0 + d)*(1.d0 + lamb)*
     $            df0dtauB(j)**3 + 
     $           80.d0*(-1.d0 + d)*(1.d0 + lamb)*df0dtauB(j)*
     $            (2.d0*d2f0dtau2B(j) + (8.d0 + 3.d0*lamb)*df0dtauB(j))*
     $            f0B(j) - 16.d0*(-1.d0 + d)*(1.d0 + lamb)*
     $            ((11.d0 + 4.d0*lamb)*d2f0dtau2B(j) + 
     $              d3f0dtau3B(j) + 
     $              (46.d0 + 34.d0*lamb + 6.d0*lamb**2)*df0dtauB(j))*
     $            f0B(j)**2 + 
     $     16.d0*(-1.d0 + d)*(1.d0 + lamb)*(2.d0 + lamb)*(3.d0 + lamb)*
     $            (4.d0 + lamb)*f0B(j)**3 + 
     $    5.d0*(-3.d0 + d)*(-2.d0 + d)**6*(-1.d0 + d*(-10.d0 + 3.d0*d))*
     $            f0B(j)**7*u1B(j)**4 + 
     $     8.d0*(-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*df0dtauB(j)*
     $            f0B(j)**4*u1B(j)*
     $            (2.d0*du1dtauB(j) + (3.d0 + lamb)*u1B(j)) - 
     $     8.d0*(-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*f0B(j)**5*
     $            (du1dtauB(j)**2 + 
     $              2.d0*(d2u1dtau2B(j) + (5.d0 + lamb)*du1dtauB(j))*
     $               u1B(j) + (9.d0 + lamb*(5.d0 + lamb))*u1B(j)**2))
     $         + 8.d0*(f0(j)*(210.d0*(-1.d0 + d)*df0dtauB(j)**3*
     $               (du1dtauB(j) + u1B(j)) - 
     $              3.d0*(-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*
     $               df0dtauB(j)*f0B(j)**4*u1B(j)**2*
     $               (du1dtauB(j) + u1B(j)) - 
     $              2.d0*(-1.d0 + d)*f0B(j)**3*
     $               (4.d0*(35.d0 + lamb*(10.d0 + lamb))*d2u1dtau2B(j) + 
     $                 6.d0*lamb*d3u1dtau3B(j) + 
     $                 4.d0*(10.d0*d3u1dtau3B(j) + d4u1dtau4B(j)) + 
     $                 (200.d0 + lamb*(80.d0 + lamb*(15.d0 + lamb)))*
     $                  du1dtauB(j) + 
     $             (6.d0 + lamb)*(16.d0 + lamb*(5.d0 + lamb))*u1B(j))
     $                - 30.d0*(-1.d0 + d)*df0dtauB(j)*f0B(j)*
     $               (4.d0*d2f0dtau2B(j)*du1dtauB(j) + 
     $                 df0dtauB(j)*
     $                  (6.d0*d2u1dtau2B(j) + 
     $                    (22.d0 + 3.d0*lamb)*du1dtauB(j)) + 
     $                 (4.d0*d2f0dtau2B(j) + 
     $                    (16.d0 + 3.d0*lamb)*df0dtauB(j))*u1B(j)) + 
     $              10.d0*(-1.d0 + d)*f0B(j)**2*
     $               (d3f0dtau3B(j)*du1dtauB(j) + 
     $                 d2f0dtau2B(j)*
     $                  (4.d0*d2u1dtau2B(j) + 
     $                    (15.d0 + 2.d0*lamb)*du1dtauB(j)) + 
     $                 2.d0*df0dtauB(j)*
     $                  ((20.d0 + 3.d0*lamb)*d2u1dtau2B(j) + 
     $                    3.d0*d3u1dtau3B(j) + 
     $                    (40.d0 + lamb*(11.d0 + lamb))*du1dtauB(j))
     $                  + ((11.d0 + 2.d0*lamb)*d2f0dtau2B(j) + 
     $                    d3f0dtau3B(j) + 
     $               2.d0*(23.d0 + lamb*(8.d0 + lamb))*df0dtauB(j))*
     $                  u1B(j)) + 
     $         (-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*f0B(j)**5*
     $               u1B(j)*(2.d0*du1dtauB(j)**2 + 
     $                 (2.d0*d2u1dtau2B(j) + 
     $                    (10.d0 + lamb)*du1dtauB(j))*u1B(j) + 
     $                 (6.d0 + lamb)*u1B(j)**2)) + 
     $           f0B(j)*(2.d0*(-1.d0 + d)*
     $               ((35.d0 + 6.d0*lamb*(5.d0 + lamb))*d2u1dtau2(j) + 
     $                 2.d0*(5.d0 + 2.d0*lamb)*d3u1dtau3(j) + 
     $                 d4u1dtau4(j) + 
     $           2.d0*(5.d0 + 2.d0*lamb)*(5.d0 + lamb*(5.d0 + lamb))*
     $                  du1dtau(j))*f0B(j)**3 - 
     $              2.d0*(-1.d0 + d)*f0B(j)**2*
     $               (4.d0*d2f0dtau2B(j)*d2u1dtau2(j) + 
     $                 4.d0*d2f0dtau2(j)*d2u1dtau2B(j) + 
     $                 40.d0*d2u1dtau2B(j)*df0dtau(j) + 
     $                 8.d0*lamb*d2u1dtau2B(j)*df0dtau(j) + 
     $                 6.d0*d3u1dtau3B(j)*df0dtau(j) + 
     $                 40.d0*d2u1dtau2(j)*df0dtauB(j) + 
     $                 18.d0*lamb*d2u1dtau2(j)*df0dtauB(j) + 
     $                 6.d0*d3u1dtau3(j)*df0dtauB(j) + 
     $                 ((15.d0 + 8.d0*lamb)*d2f0dtau2B(j) + 
     $                    d3f0dtau3B(j) + 
     $                    2.d0*(40.d0 + lamb*(40.d0 + 9.d0*lamb))*
     $                     df0dtauB(j))*du1dtau(j) + 
     $                 15.d0*d2f0dtau2(j)*du1dtauB(j) + 
     $                 3.d0*lamb*d2f0dtau2(j)*du1dtauB(j) + 
     $                 d3f0dtau3(j)*du1dtauB(j) + 
     $                 80.d0*df0dtau(j)*du1dtauB(j) + 
     $                 30.d0*lamb*df0dtau(j)*du1dtauB(j) + 
     $                 3.d0*lamb**2*df0dtau(j)*du1dtauB(j) + 
     $                 ((11.d0 + 3.d0*lamb)*d2f0dtau2(j) + 
     $                    d3f0dtau3(j) + 
     $              (46.d0 + lamb*(22.d0 + 3.d0*lamb))*df0dtau(j))*
     $                  u1B(j)) + 
     $              10.d0*(-1.d0 + d)*f0B(j)*
     $               (6.d0*d2u1dtau2B(j)*df0dtau(j)*df0dtauB(j) + 
     $                 3.d0*d2u1dtau2(j)*df0dtauB(j)**2 + 
     $                 df0dtauB(j)*
     $                  (2.d0*d2f0dtau2B(j) + 
     $                    (11.d0 + 6.d0*lamb)*df0dtauB(j))*du1dtau(j)
     $                  + 2.d0*d2f0dtau2B(j)*df0dtau(j)*
     $                  du1dtauB(j) + 
     $                 2.d0*d2f0dtau2(j)*df0dtauB(j)*du1dtauB(j) + 
     $                 22.d0*df0dtau(j)*df0dtauB(j)*du1dtauB(j) + 
     $                 4.d0*lamb*df0dtau(j)*df0dtauB(j)*
     $                  du1dtauB(j) + 
     $                 2.d0*(d2f0dtau2(j)*df0dtauB(j) + 
     $                    df0dtau(j)*
     $                     (d2f0dtau2B(j) + 
     $                       2.d0*(4.d0 + lamb)*df0dtauB(j)))*u1B(j))
     $               - (-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*
     $               f0B(j)**5*u1B(j)*
     $               (2.d0*du1dtau(j)*du1dtauB(j) + 
     $                 (d2u1dtau2(j) + (5.d0 + 2.d0*lamb)*du1dtau(j))*
     $                  u1B(j)) + 
     $        (-3.d0 + d)*(-2.d0 + d)**3*(-1.d0 + 4.d0*d)*f0B(j)**4*
     $               u1B(j)**2*
     $               (df0dtauB(j)*du1dtau(j) + 
     $                 df0dtau(j)*(du1dtauB(j) + u1B(j))) - 
     $              30.d0*(-1.d0 + d)*df0dtauB(j)**2*
     $               (df0dtauB(j)*du1dtau(j) + 
     $                 3.d0*df0dtau(j)*(du1dtauB(j) + u1B(j))))))/
     $      (128.d0*(-1.d0 + d)**2*(1.d0 + d)*f0B(j)**8)

         v5(j) = u5(j)
      end do


C     Put together the expansion around x=0. 
      do j=1,ny
         f(j) = f0(j)  + x**2 * f2(j) 
     $    + x**4 * f4(j)
         u(j) = x * u1(j) + x**2 * u2(j) + x**3 * u3(j) 
     $           + x**4 * u4(j) + x**5 * u5(j)
         v(j) = x * v1(j) + x**2 * v2(j) + x**3 * v3(j) 
     $           + x**4 * v4(j) + x**5 * v5(j)
         ia2(j) = x**2 *ia22(j) 
     $    + x**4 * ia24(j)
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
            write(11,*) f0(j), f2(j), f4(j)
         end do
      close(11)
      open(unit=11,file='uleftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) u1(j), u2(j), u3(j), u4(j), u5(j)
         end do
      close(11)
      open(unit=11,file='vleftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) v1(j), v2(j), v3(j), v4(j), v5(j)
         end do
      close(11)
      open(unit=11,file='ia2leftcoeffs.junk',status='new')
         do j=1,ny
            write(11,*) ia22(j), ia24(j)
         end do
      close(11)
      
      end if

      call yfromfields(ny, u, v, f, y)
    

      return
      end
