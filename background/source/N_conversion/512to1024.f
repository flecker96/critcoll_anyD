C***********************************************************************
C
C     This program reads in .out files from the converged result at N=128
C     and pads them with zero modes such that the result can be used as
C     initial data for an N=512 computation
C
C***********************************************************************

      program double_dat

C     *******************************
C     **** Variable declarations ****
C     *******************************

      implicit none
      
C     tau related variables
      integer ny, nymax, j
      include '../nymax.inc'

      double precision Delta, fc(nymax), Up(nymax), psic(nymax),
     $       fcF(nymax), UpF(nymax), psicF(nymax),
     $       fout(nymax), uout(nymax), psiout(nymax)
     

      ny = 1024

C     *******************
C     **** Free data ****
C     *******************

C     Read in free functions from dat files. 
      open(unit=10, file='512/Delta.out', status='old')
      read(10,*) Delta
      close(10)
      open(unit=10, file='512/fc.out', status='old')
      do j=1,ny
         read(10,*) fc(j)
      end do
      close(10)
      open(unit=10, file='512/psic.out', status='old')
      do j=1,ny
         read(10,*) psic(j)
      end do
      close(10)    
      open(unit=10, file='512/Up.out', status='old')
      do j=1,ny
         read(10,*) Up(j)
      end do
      close(10)
      
      write(6,*) 'INFO: Upsampling initial data.'
      write(6,'(A13,I4,A9,I4)') ' INFO: Ntin = ',ny/2,', Nout = ',ny

      do j=1,ny
         fcF(2*j-1) = fc(j) / ny
         psicF(2*j-1) = psic(j) / ny
         UpF(2*j-1) = Up(j) / ny
         fcF(2*j) = 0.d0
         psicF(2*j) = 0.d0
         UpF(2*j) = 0.d0
      end do

      call four1(fcF, ny, 1)
      call four1(psicF, ny, 1)
      call four1(UpF, ny, 1)

C     Double modes once
      call double(fcF, fcF, ny)
      call double(psicF, psicF, ny)
      call double(UpF, UpF, ny)

      call four1(fcF, 2*ny, -1)
      call four1(psicF, 2*ny, -1)
      call four1(UpF, 2*ny, -1)
      
C     take only real parts for functions      
      do j=1,2*ny
         fout(j) = fcF(2*j-1)
         psiout(j) = psicF(2*j-1)
         uout(j) = UpF(2*j-1)
      end do


C     **********************
C     *** Write to file ****
C     **********************

      call system('mkdir 1024')

      open(unit=10, file='1024/Delta.dat', status='new')
      write(10,99) Delta
      close(10)

      open(unit=10, file='1024/fc.dat', status='new')
      do j=1,2*ny
         write(10,99) fout(j)
      end do
      close(10)

      open(unit=10, file='1024/psic.dat', status='new')
      do j=1,2*ny
         write(10,99) psiout(j)
      end do
      close(10)     

      open(unit=10, file='1024/Up.dat', status='new')
      do j=1,2*ny
         write(10,99) uout(j)
      end do
      close(10)

 99   format(G24.16)
      end



      
