C     nu-stage Implicit Runge-Kutta (2nu-order), solved by iteration

      subroutine implicit_steplin(ny, d, Delta, lamb, 
     $        xin, yin, xout, yout,
     $        nu, crit, maxits, itsreach, 
     $        uBcoll, vBcoll, fBcoll, ia2Bcoll)

      implicit none
      integer ny, nu, maxits, itsreach
      double precision xin, yin(ny), xout, yout(ny), crit, d, lamb
      external derivslin

      integer nymax, numax
      include '../nymax.inc'

      parameter(numax=3)
      double precision xi(nymax,numax), xiold(nymax,numax), 
     $    uBcoll(nymax,numax), vBcoll(nymax,numax), 
     $    fBcoll(nymax,numax), ia2Bcoll(nymax,numax),
     $    norminfdiff, norm2diff, pointdiff,
     $    crit10
      double precision dx, Delta, x(numax), f(nymax,numax)
      double precision a(numax,numax), b(numax), c(numax)
      double precision ef(nymax), u(nymax), v(nymax), ia2(nymax)
      integer i, ii, j, its

      double precision y(nymax)

      logical paus
      
      paus = .FALSE.

      if (ny .gt. nymax) stop 'nymax too small in implicit_step'

C     Define IRK nu=1,2,3
      if (nu.eq.1) then
         a(1,1) = 0.5d0
         b(1) = 1.d0
         c(1) = 0.5d0
      else if (nu.eq.2) then
         a(1,1) = 0.25d0
         a(1,2) = 0.25d0 - 0.5d0 / sqrt(3.d0)
         a(2,1) = 0.25d0 + 0.5d0 / sqrt(3.d0)
         a(2,2) = 0.25d0
         b(1) = 0.5d0
         b(2) = 0.5d0
         c(1) = 0.5d0 - 0.5d0 / sqrt(3.d0)
         c(2) = 0.5d0 + 0.5d0 / sqrt(3.d0)
      else if (nu.eq.3) then
         a(1,1) = 5.d0/36.d0
         a(1,2) = 2.d0/9.d0 - 1.d0/sqrt(15.d0)
         a(1,3) = 5.d0/36.d0 - 0.5d0/sqrt(15.d0)
         a(2,1) = 5.d0/36.d0 + sqrt(15.d0)/24.d0
         a(2,2) = 2.d0/9.d0
         a(2,3) = 5.d0/36.d0 - sqrt(15.d0)/24.d0
         a(3,1) = 5.d0/36.d0 + 0.5d0/sqrt(15.d0)
         a(3,2) = 2.d0/9.d0 + 1.d0/sqrt(15.d0)
         a(3,3) = 5.d0/36.d0
         b(1) = 5.d0/18.d0
         b(2) = 4.d0/9.d0
         b(3) = 5.d0/18.d0
         c(1) = 0.5d0 - sqrt(15.d0)/10.d0
         c(2) = 0.5d0
         c(3) = 0.5d0 + sqrt(15.d0)/10.d0
      else
         write(*,*) 'IRK nu=', nu, ' not implemented in implicit_step.'
         stop
      end if
      
C     Start new calculation

      do j=1,ny
         y(j) = yin(j)
      end do


C     Collocation points
      dx = ( xout - xin ) 
      do i=1,nu
         x(i) = xin + dx * c(i)
      end do

     
C     Zeroth-order estimates
      do j=1,ny
         do i=1,nu
            xi(j,i) = y(j)
         end do
      end do
      

      do its=1,maxits
C        Store values of stages from previous iteration
         do j=1,ny
            do i=1,nu
               xiold(j,i) = xi(j,i)
            end do
         end do
C        Evaluate derivatives f at stages xi
         do i=1,nu
            call derivslin(ny, d, Delta, lamb, xi(1,i), x(i), f(1,i),
     $        uBcoll(1,i), vBcoll(1,i), fBcoll(1,i), ia2Bcoll(1,i))
C          write(6,*) xi(i,1)
         end do         
                   
C        Calculate new stages
         do j=1,ny
            do i=1,nu
               xi(j,i) = y(j)
               do ii=1,nu
                  xi(j,i) = xi(j,i) + dx * a(i,ii) * f(j,ii)
               end do
            end do
         end do
C        2-Norm of difference
         norm2diff = 0.d0
         do j=1,ny
            do i=1,nu
               norm2diff = norm2diff + (xiold(j,i)-xi(j,i))**2
            end do
         end do
         
         norm2diff = sqrt( norm2diff / nu / ny )
         
C        Decision point
         if (its .le. maxits/2) then
            crit10 = crit
         else
            crit10 = crit * 10.d0 * (2.d0*dble(its)/dble(maxits) - 1.d0)
         end if 
         
         if (norm2diff .lt. crit10 ) then
C           Check inf-norm of difference
C           By definition we have that inf-norm <= sqrt(ny) * 2-norm,
C           but for large ny this could be too large a difference.
C           Therefore we inforce a maximum of a factor 10.
            norminfdiff = 0.d0
             
            do j=1,ny
               do i=1,nu
                 pointdiff = abs(xiold(j,i)-xi(j,i))
                 norminfdiff = max(norminfdiff, pointdiff)
               end do
            end do
          
            if (norminfdiff .lt. 10.d0*crit10 ) goto 2
            
         end if
      end do
      write(6,*) 'implicitstep does not converge at x= ', xin 
C      stop 
      paus = .TRUE.
C     The program only branches here when the iteration has converged.
 2    continue
      
      itsreach = max( its, itsreach )
      
C     Evaluate new value of y
      do j=1,ny
         do i=1,nu
            y(j) = y(j) + dx * b(i) * f(j,i)
         end do
         yout(j) = y(j)
      end do

      if (paus) then
      
      call fieldsfromy(ny, d, Delta, lamb, yout, xout,
     $           u, v, ef, ia2)
      
      open(unit=10,file='ytest.junk',status='new')
      
         do i=1,ny
            write(10,*) ef(i), u(i), v(i)
         end do
      close(10)
      stop
      end if


      return

      end

