
      subroutine myunpack(n3, vec, ny, odd1, odd2, even)

      implicit none
      integer n3, ny
      double precision vec(n3), odd1(ny), odd2(ny), even(ny)

      integer nymax, j
      include '../nymax.inc'
      double precision odd1F(2*nymax), odd2F(2*nymax), evenF(2*nymax)

      do j=1,2*ny
         odd1F(j) = 0.d0
         odd2F(j) = 0.d0
         evenF(j) = 0.d0
      end do
      do j=1,n3/6
         odd1F(4*j-1) = vec(2*j-1       )
         odd1F(4*j  ) = vec(2*j         ) 
         odd1F(ny-4*j+3) = odd1F(4*j-1)
         odd1F(ny-4*j+4) = - odd1F(4*j)
         odd2F(4*j-1) = vec(2*j-1+  n3/3)
         odd2F(4*j  ) = vec(2*j  +  n3/3)
         odd2F(ny-4*j+3) = odd2F(4*j-1)
         odd2F(ny-4*j+4) = - odd2F(4*j)
         evenF(4*j-3) = vec(2*j-1+2*n3/3)
         evenF(4*j-2) = vec(2*j  +2*n3/3)
         evenF(ny-4*j+5) = evenF(4*j-3)
         evenF(ny-4*j+6) = - evenF(4*j-2)
      end do

C     Set the component to zero again where the high freq cos was stored      
      evenF(2)= 0.d0
      
      call double(odd1F, odd1F, ny/2)
      call double(odd2F, odd2F, ny/2)
      call double(evenF, evenF, ny/2)
      
C     High frequency cosine
      evenF(2*n3/3 + 1) = vec(2*n3/3+2)/2.d0
      evenF(2*n3   + 1) = vec(2*n3/3+2)/2.d0
      
      
      call four1(odd1F, ny, -1)
      call four1(odd2F, ny, -1)
      call four1(evenF, ny, -1)
      
C     take only reals parts for functions      
      do j=1,ny
         odd1(j) = odd1F(2*j-1)
         odd2(j) = odd2F(2*j-1)
         even(j) = evenF(2*j-1)
      end do


      return
      end
