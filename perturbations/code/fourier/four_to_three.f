      subroutine four_to_three(n4, y4, n3, y3)

      implicit none
      integer n4, n3
      double precision y4(n4), y3(n3)

      integer i

      if (n4 .ne. n3*4/3) stop 'n4 ne 4/3 n3 in four_to_three'

C     Odd frequencies.
      do i = 0, n4 / 4 - 1
         y3(2 * i + 1) = y4(4 * i + 3)
         y3(2 * i + 2) = y4(4 * i + 4)
      end do

C     Imaginary part of zero frequency. (Compare i = 0 below.)
      y3(n4 / 2 + 1) = y4(2)

C     Imaginary part of high frequency cosine. (Compare i = n4/8 below.)
      y3(n4 / 2 + 2) = y4(n4 / 2 + 2)

C     Imaginary part of even frequencies: sine and cosine.
      do i = 1, n4 / 8 - 1
         y3(n4 / 2 + 2 * i + 1) = 
     $        0.5d0 * (y4(4 * i + 1) - y4(n4 - 4 * i + 1))
         y3(n4 / 2 + 2 * i + 2) = 
     $        0.5d0 * (y4(4 * i + 2) + y4(n4 - 4 * i + 2))
      end do

C    Structure of array y3: 
C    
C    ( Re(u_1)-Im(v_1),Im(u_1)-Re(v_1), Re(u_3)-Im(v_3),Im(u_3)-Re(v_3),
C      ... , Im(u_(n4/2 - 1))-Re(v_(n4/2 - 1)), Re(f_0), Re(f_n4/4) ,
C     -Im(f_2), Re(f_2), -Im(f_4), Re(f_4), 
C      ... ,-Im(f_(n4/4-2)), Re(f_(n4/4 - 2))


      return
      end


