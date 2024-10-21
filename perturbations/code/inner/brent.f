      subroutine brent(func,x1,x2,tol)

      integer ITMAX
      double precision zbrent,tol,x1,x2,func,EPS
      external func
      parameter (ITMAX=100,EPS=3.d-8)
      integer iter
      double precision a,b,c,d,e,fa,fb,fc,p,q,r,
     $     s,tol1,xm

      a=x1
      b=x2
      fa=func(a)
      fb=func(b)

      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
        stop 'Root not bracketed in brent.'
      end if

      c=b
      fc=fb

      do iter=1,ITMAX

C       If c is on same side as b, set c=a       
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.)) then
            c=a
            fc=fa
            d=b-a
            e=d
        endif

C       Make sure that b is the better guess, switch with a if not the case
        if(abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        endif

C       xm is (signed) distance between b and midpoint
        tol1=2.d0*EPS*abs(b) + tol/2.d0
        xm=(c-b) / 2.d0

        if(abs(xm).le.tol1 .or. fb.eq.0.) then
            zbrent=b
            goto 1
        endif

        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
            s=fb/fa
            if(a.eq.c) then
                p=2.d0*xm*s
                q=1.d0 - s
            else
                q=fa/fc
                r=fb/fc
                p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
                q=(q-1.d0)*(r-1.d0)*(s-1.d0)
            endif

C           Check whether in bounds.
            if(p.gt.0.) then 
                q=-q
            end if

            p=abs(p)
            
            if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
C               Accept interpolation.
                e=d
                d=p/q
            else
C               Interpolation failed, use bisection.
                d=xm
                e=d
            endif
        else
C       Bounds decreasing too slowly, use bisection.
            d=xm
            e=d
        endif

C       Move last best guess to a.
        a=b
        fa=fb

C       Evaluate new trial root.
        if(abs(d) .gt. tol1) then
            b=b+d
        else
            b=b+sign(tol1,xm)
        endif

        fb=func(b)

      end do

      stop 'brent exceeding maximum iterations'
      
 1    zbrent=b
      return
      END