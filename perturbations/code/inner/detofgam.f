      subroutine detofgam(ny, nx, d, gam, ileft, iright, imid,
     $      n3, in0, in0B, outevery, xxp, prec_irk, debug, 
     $      uB, vB, fB, ia2B, ivarread, its, det, 
     $      doutdin, eps, fac, newfac)
     
      implicit none
      
C     tau related variables
      integer ny, n3, nymax, n3max, j
      include '../nymax.inc'
      parameter(n3max = 3*nymax/4)
      double precision in0(n3max), out0(n3max),
     $      lamb, gam
 
C     x related variables
      integer nx, ileft, imid, iright, outevery,
     $     nxmax, itsreach
      parameter(nxmax = 300001)
      double precision prec_irk, xxp(nxmax)
     
      logical debug, newfac

C     Background variables
      double precision  
     $     uB(nymax, nxmax), vB(nymax,nxmax), 
     $     fB(nymax,nxmax), ia2B(nymax, nxmax),
     $     in0B(n3max)

C     Newton method variables
      integer its, indx(n3max), ivar, maxits, ivarread
      parameter(maxits=100)
      double precision in1(n3max), out1(n3max),
     $     doutdin(n3max,n3max),
     $     eps, d, err, fac,
     $     det,
     $     log_eig_sum
         
      lamb = 1.d0 / gam

C     Base shot
      itsreach = 1

      call shootlin(ny, nx, d, lamb, ileft, iright, imid, 0,
     $      n3, in0, in0B, out0, outevery, xxp, 2, 
     $      prec_irk, debug, itsreach, 
     $      uB, vB, fB, ia2B)

C     Calculate mismatch (just for comparison, not important here)
      err = 0.d0
      do j=1,n3
      err = err + out0(j)**2
      end do
      err = sqrt(err/n3)
      write(*,*) 'Perturbation mismatch: ', err
      
C     Compute Jacobian
      do ivar=ivarread+1,n3

C           Store input free data
            do j=1,n3
                  in1(j) = in0(j)
            end do
C           Change variable ivar
            in1(ivar) = in1(ivar) + eps

C           Output info
            write(6,*) '  Max its reached:', itsreach, '. Varying', ivar
            call flush(6)

C           Shoot new free data
            itsreach = 1
            call shootlin(ny, nx, d, lamb, ileft, iright, imid, ivar,
     $       n3, in1, in0B, out1, 0, xxp, 2, prec_irk, debug,  
     $       itsreach, uB, vB, fB, ia2B)

C           Calculate derivative with respect to variable ivar
            do j=1,n3
                  doutdin(j,ivar) = (out1(j) - out0(j)) / eps
            end do
            
            call output_step(n3, its+1, ivar, doutdin(1,ivar))

      end do
      

C     Here, the determinant is computed. To avoid overflow the entries of the 
C     matrix are rescaled by 10^(-log(det(A))/n3)
C     LU decomposition of doutdin, stored in the same space
      call ludcmp(doutdin, n3, n3max, indx, det)
      
      if (newfac) then
            log_eig_sum = 0.d0

            do j = 1, n3
                  log_eig_sum = log_eig_sum + log10(abs(doutdin(j,j)))
            end do

            fac = 10.d0**( - (log_eig_sum) / dble(n3))
      end if

      do j = 1, n3
            det = det * doutdin(j,j) * fac
      end do


      return
      end