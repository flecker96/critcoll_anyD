
      subroutine mypack(n3, vec, ny, odd1, odd2, even)

      implicit none
      integer n3, ny
      double precision vec(n3), odd1(ny), odd2(ny), even(ny)

      integer nymax, j
      include '../nymax.inc'
      double precision odd1F(2*nymax), odd2F(2*nymax), evenF(2*nymax)

      do j=1,ny
         odd1F(2*j-1) = odd1(j) / ny
         odd2F(2*j-1) = odd2(j) / ny
         evenF(2*j-1) = even(j) / ny
         odd1F(2*j) = 0.d0
         odd2F(2*j) = 0.d0
         evenF(2*j) = 0.d0
      end do
      call four1(odd1F, ny, 1)
      call four1(odd2F, ny, 1)
      call four1(evenF, ny, 1)
      
     
C     Eliminate upper half of frequencies (anti-aliasing). 
      call halve(odd1F, odd1F, ny/2)
      call halve(odd2F, odd2F, ny/2)
      call halve(evenF, evenF, ny/2)
     
      
C     Use the symmetries of the three functions to decrease
C     the number of real variables to n3 = 3*ny/4
C     (Recall that ny is already the doubled number of modes)
      do j=1,n3/6
         vec(2*j-1       ) = odd1F(4*j-1)
         vec(2*j         ) = odd1F(4*j  )
         vec(2*j-1 +  n3/3) = odd2F(4*j-1)
         vec(2*j   +  n3/3) = odd2F(4*j  )
         vec(2*j-1 + 2*n3/3) = evenF(4*j-3)
         vec(2*j   + 2*n3/3) = evenF(4*j-2)
      end do
      
C     High frequency cosine replaces the imaginary part of 0 frequency
C     of the even function, which vanishes.

      vec(2*n3/3+2)=evenF(2*n3/3+1)
      
      return
      end
