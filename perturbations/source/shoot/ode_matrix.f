C     Compute the coefficient matrix for the linearized ODE - system 
C     in Fourier space

      subroutine jacobian(n3, ny, d, lamb, Delta, x, y3, jac,
     $     uB, vB, fB, ia2B)

      implicit none
      integer ny, n3
      double precision x, y3(n3), lamb, Delta, d, eps
      external derivslin

      integer nymax, numax
      include '../nymax.inc'

     
      double precision jac(n3,n3), uB(nymax), vB(nymax),
     $                 fB(nymax), ia2B(nymax), test(nymax)
      double precision y3new(nymax), dy3old(nymax), dy3new(nymax)

      integer i, ii, j, its

      if (ny .gt. nymax) stop 'nymax too small in implicit_step'

      call derivslin(n3, ny, d, Delta, lamb, y3, x, dy3old,
     $              uB, vB, fB, ia2B) 

      do i = 1,ny
         y3new(i) = y3(i)
      end do

      do j = 1,n3
         y3new(j) = y3(j) + eps
         call derivslin(n3, ny, d, Delta, lamb, ynew, x, dynew,
     $              uB, vB, fB, ia2B) 
         do i = 1,ny
            mevo(i,j) = (dynew(i) - dyold(i)) / eps
         end do
         ynew(j) = y(j)      
      end do
          
      
          
C      do i = 1,ny
C         test(i) = 0.d0
C         do j = 1,ny
C            test(i) = test(i) + mevo(i,j) * y(j)
C         end do
C      
C         write(6,*) dyold(i) - test(i)
C      end do

      return

      end

