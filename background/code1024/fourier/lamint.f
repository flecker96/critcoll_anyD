C     From f(x) = sum f_k exp(-ikx) calculates the generalised principal
C     function F(lambda,x) = sum f_k / (lam - i k) exp(-ikx). The 
c     special case lambda = 0, F(0,x) = sum_(k.ne.0) f_k / (-ik) exp(ikx).
C     Also returns the constant part f0 of f.

      subroutine lamint(n, f, period, lam, f0, bigf)

      integer n, nymax
      include '../nymax.inc'
      double precision f(n), period, lam, f0, bigf(n),
     $     caux(2*nymax)

      integer i
      double precision k1, retemp, imtemp

      if (n .gt. nymax) stop 'nymax too small in lamint'
C     k1 is the fundamental frequency
      k1 = 6.28318530717959d0 / period

C     From f(x) calculate f0 and bigf(x), by dividing f_k by (lam-ik).
      do i=1,n
         caux(2*i-1) = f(i) / dble(n)
         caux(2*i) = 0.d0
      end do
      call four1(caux, n, 1)

      f0 = caux(1) 

C     Zero frequency
      if (lam .eq. 0.d0) then
         caux(1) = 0.d0
         caux(2) = 0.d0
      else
         caux(1) = caux(1) / lam
         caux(2) = caux(2) / lam
      end if

C     Positive frequencies:
      do i=2,n/2
         retemp = lam * caux(2*i-1) - dble(i-1) * k1 * caux(2*i)
         imtemp = lam * caux(2*i) + dble(i-1) * k1 * caux(2*i-1)
         caux(2*i-1) = retemp / (lam**2 + dble(i-1)**2 * k1**2)
         caux(2*i) = imtemp / (lam**2 + dble(i-1)**2 * k1**2)
      end do

C     High frequency treated as a cosine:
      caux(n+1) = caux(n+1) * lam / (lam**2 + dble(n/2)**2 * k1**2)
      caux(n+2) = caux(n+2) * lam / (lam**2 + dble(n/2)**2 * k1**2)

C     Negative frequencies.
      do i=n/2+2,n
         retemp = lam * caux(2*i-1) - dble(i-1-n) * k1 * caux(2*i)
         imtemp = lam * caux(2*i) + dble(i-1-n) * k1 * caux(2*i-1)
         caux(2*i-1) = retemp / (lam**2 + dble(i-1-n)**2 * k1**2)
         caux(2*i) = imtemp / (lam**2 + dble(i-1-n)**2 * k1**2)
      end do

C     Fourier transform back.
      call four1(caux,n,-1)
      do i=1,n
         bigf(i) = caux(2*i-1)
      end do

      return
      end








