      subroutine three_to_four(n3, y3, n4, y4)

      implicit none
      integer n3, n4
      double precision y3(n3), y4(n4)

      integer i

      if (n4 .ne. n3*4/3) stop 'n4 ne 4/3 n3 in three_to_four'

C     Odd frequencies.
      do i = 0, n4/4 - 1
         y4(4 * i + 3) = y3(2 * i + 1)
         y4(4 * i + 4) = y3(2 * i + 2)
      end do

C     Zero frequency purely imaginary.
      y4(1) = 0.d0
      y4(2) = y3(n4 / 2 + 1)

C     High frequency cosine purely imaginary.
      y4(n4 / 2 + 1) = 0.d0
      y4(n4 / 2 + 2) = y3(n4 / 2 + 2)

C     Even frequencies purely imaginary.
      do i = 1, n4 / 8 - 1
         y4(4 * i + 1) = y3(n4 / 2 + 2 * i + 1)
         y4(n4 - 4 * i + 1) = - y3(n4 / 2 + 2 * i + 1)
         y4(4 * i + 2) = y3(n4 / 2 + 2 * i + 2)
         y4(n4 - 4 * i + 2) = y3(n4 / 2 + 2 * i + 2)
      end do

      return
      end

C     Discrete Fourier transform
C                                         n k
C     f_n = Sum(k=1...N-1) F_k exp(2 Pi i ---)
C                                          N
C     Therefore 
C     
C     Re f_n = 1/2 Sum [F_k + F_(N-k)*] exp...
C     
C     This must vanish for even k, and so for even k
C
C     F_k     =  a_k + i b_k
C    
C     F_(N-k) = -a_k + i b_k
C
C     In y3, we have 
C     Re F_1, Im F_1, Re F_3, Im F_3, ..., Re F_N-1, Im F_N-1 
C     followed by  
C     b_0, a_2, b_2, ..., a_N/2-2, b_N/2-2, b_N/2
C     while
C     a_0 = 0 and a_N/2 = 0, from the condition above.
