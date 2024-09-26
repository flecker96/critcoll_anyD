      SUBROUTINE FOUR1(data,nn,isign)
C--------------------------------------------------------------------------
      IMPLICIT NONE
      integer nn,isign
      double precision data(2*nn)
C--------------------------------------------------------------------------
C     From Press et al, Fourier transform of complex function in nn points.
C     nn must be a power of 2. isign=-1 for inverse Fourier transform. 
C     (Inverse transform must be divided by nn in calling routine.)
C--------------------------------------------------------------------------
      integer n,j,i,m,mmax,istep
      double precision wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
C--------------------------------------------------------------------------

      n=2*nn
      j=1
      do i=1,n,2
         if (j.gt.i) then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         end if
         m=n/2
         do while ((m.ge.2).and.(j.gt.m)) 
            j=j-m
            m=m/2
         end do
         j=j+m
      end do
      mmax=2
      do while (n.gt.mmax) 
         istep=2*mmax
         theta=6.28318530717959d0/dble(isign*mmax)
         wpr=-2.d0*dsin(0.5d0*theta)**2
         wpi=dsin(theta)
         wr=1.d0
         wi=0.d0
         do m=1,mmax,2
            do i=m,n,istep
               j=i+mmax
               tempr=wr*data(j)-wi*data(j+1)
               tempi=wr*data(j+1)+wi*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         end do
         mmax=istep
      end do
      
      return
      end
              

