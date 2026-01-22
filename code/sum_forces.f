c      This program, QFILTER modifies classical trajectories according 
c      to the GSTA method:
c      Dénes Berta, Dávid Ferenc, Imre Bakó and Ádám Madarász 
c      "Nuclear Quantum Effects from the Analysis of Smoothed Trajectories: 
c      Pilot Study for Water"
c      https://doi.org/10.1021/acs.jctc.9b00703
c      email: madarasz.adam@ttk.hu
c     
c      QFILTER works interactvely, just start the binary file, then
c      1. you need to give the name of the xyz trajectory file
c      2. the program asks for the timestep between two snapshots
c      3. the temperature of the simulation must be entered
c      

      program qfilter

c These parameters determines the numerical integrals:
      DOUBLE PRECISION ftstep,kercut,precint

      parameter (ftstep=1.0d-4)    ! step size in the Fourier transform
      parameter (kercut=1.0d-6)    ! cutoff value of the kernel function
      parameter (precint=1.0d-11)  ! precision of the integral

      DOUBLE PRECISION avogadro,boltzmann,planck

      parameter (avogadro=6.02214129d+23)
      parameter (boltzmann=0.831446215d0)
      parameter (planck=6.62606957d-34)
      integer gnum        ! number of table for the weight function
      integer ndec,ablak,nblock,natoms
      integer i,j,k,m,ipos,tav,midpos,nst,glim,nstmax
      integer ierror,iline
      integer pbead            ! P: number of replacs or beads
      integer nspacx,nspacy,nspacz
      DOUBLE PRECISION pbcx,pbcy,pbcz
      DOUBLE PRECISION hpbcx,hpbcy,hpbcz
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: bx(:),by(:),bz(:)
      DOUBLE PRECISION, allocatable :: cx(:),cy(:),cz(:)
      character*120 bxyzfile,cxyzfile,outxyzfile,stratoms
      character(300), allocatable :: title(:)  
      character(5), allocatable :: nam(:)  

      DOUBLE PRECISION, allocatable :: gx(:,:),gy(:),hy(:)
      DOUBLE PRECISION, allocatable :: hapy(:)
      DOUBLE PRECISION, allocatable :: xkxk(:),gynew(:),gyold(:)

      DOUBLE PRECISION xa,xb,ya,yb,ipgy
      DOUBLE PRECISION gmax,dx,loc,poc,relx,egyik,masik

      DOUBLE PRECISION absx,sumu,sumb,oseg,candy,step

      DOUBLE PRECISION hiba,zaj,ero

      common pbead

      DOUBLE PRECISION gweight

      DOUBLE PRECISION outdx,outgmax

      LOGICAL :: gx_file_exist,kernel_file_exist

      character*120 gx_file_name,kernel_file_name

      INTEGER outgnum

      kernel_file_name="kernel.dat"

      PI=4.D0*DATAN(1.D0)

c
c     get the base name of user specified input structures
c
      write (*,10)
  10    format (/,' Name of the first xyz file: ',$)
      read (*,20)  bxyzfile
   20    format (a120)

      write (*,*) bxyzfile

      write (*,12)
  12    format (/,' Name of the second xyz file: ',$)
      read (*,20)  cxyzfile

      write (*,*) cxyzfile

      write (*,13)
  13    format (/,' Name of the output xyz file: ',$)
      read (*,20)  outxyzfile

      write (*,*) outxyzfile


c
c     get the number of frames
c

      OPEN (50, file = bxyzfile)
      read(50,'(i10)',iostat=ierror) natoms

      if ( ierror .ne. 0 ) then

            write (*,*) "ERROR: Incorrect xyz file format. ",
     +          "The first line is not an integer ",
     +          "that should mean the number of atoms. Program stops."
            stop

      endif

      i = 1
      DO
       READ (50,*, END=11)
       i = i + 1
      END DO
   11 CLOSE (50)

      if (mod(i,(natoms+2)) .ne. 0) then

            write (*,*) "WARNING: The number of lines does not ",
     +          "match the number of atoms."
            stop

      endif 

      nblock=i/(natoms+2)

      write(*,*) 'Number of atoms: ',natoms
      write(*,*) 'Number of input frames: ',nblock

c      write (*,80)
c   80 format (/,' Length of the simulation box ',
c     &              ' in Angstrom:  ',$)
c   90 format (f20.0)
c      read (*,90) pbcx
c      write (*,*) pbcx

c      pbcy=pbcx
c      pbcz=pbcx


c      hpbcx=pbcx/2
c      hpbcy=pbcy/2
c      hpbcz=pbcz/2

c      write (*,81)
c   81 format (/,' Temperature for filtration',
c     &             ' in Kelvin:  ',$)
c   91    format (f20.0)
c         read (*,91) temp
c         write (*,*) temp
c      write(*,*) 'Calculation of the kernel function:'

c      deltat=2*PI/(avogadro*planck/boltzmann/temp/timestep*1E11)

c Reading the kernel function

c      INQUIRE(FILE=kernel_file_name, EXIST=kernel_file_exist)

c      if (kernel_file_exist) then

c        OPEN (51, file = kernel_file_name)

c         j = 0
c         DO
c           READ (51,*, END=14) kernel(j)
c           j = j + 1
c         END DO
c   14 CLOSE (51)

c      else

c          write(*,*) "kernel.dat file is not found. Program stops."
c          stop

c      endif

      kernel(0)=1.0d0

      ndec=0

      ablak=2*ndec+1

c Check the integral of the kernel function

c      sum=kernel(0)/2.0d0

c      do j=1,ndec

c          sum=sum+kernel(j)

c      end do

c      sum=sum*2.0d0

c      write(*,*) "The integral of the kernel function is: ",sum
c      write(*,*) "before it is normalized to 1."
c Normalize the kernel function to 1.0

c      do j=0,ndec

c          write(*,*) kernel(j)

c          kernel(j)=kernel(j)/sum

c      end do

      write(*,*) 'Starting filtration'
c     End of the calculation of the kernel function



c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock))
      allocate (nam(natoms))
      allocate (bx(natoms))
      allocate (by(natoms))
      allocate (bz(natoms))

      allocate (cx(natoms))
      allocate (cy(natoms))
      allocate (cz(natoms))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      open (unit=50,file=bxyzfile,status='old',action='read')  

      open (unit=55,file=cxyzfile,status='old',action='read')  

      open (unit=60,file=outxyzfile,status='unknown',action='write')  

      write(*,*) 'Name of the outputfile: ', outxyzfile
      write(*,*) 'Number of filtered frames: ', nblock-ablak+1

      iline=0

c
c     cycle over all pairs of snapshot frame blocks
c
      do i = 1, nblock

         read(50,'(a)') stratoms
         read(50,15) title(i)
         read(55,'(a)') stratoms
         read(55,15) title(i)
   15    format (a300)

         iline=iline+2

         do j=1,natoms
            iline=iline+1
            read(50,*,iostat=ierror) nam(j),bx(j),by(j),bz(j)  
            read(55,*,iostat=ierror) nam(j),cx(j),cy(j),cz(j)  

            if ( ierror .ne. 0 ) then

                write (*,*) "ERROR: Problem reading files "
                write (*,*) "in line ",iline
                write (*,*)  "Program stops."
            stop

            endif

         end do  
 
         ipos=mod(i-1,ablak)+1
         midpos=mod(i-1-ndec,ablak)+1

c
c     write output
c

      write(*,*) 'Writing frame', i-2*ndec
      
            write(60,'(a)') trim(stratoms)

            write(60,'(a)') trim(title(i-ndec))

            do k=1,natoms
               write(60,25) nam(k),bx(k)+cx(k),by(k)+cy(k),bz(k)+cz(k)
   25    format (a3,3f16.10)
            end do  

      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (bx)
      deallocate (by)
      deallocate (bz)

      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)


      close(50)
      close(55)
      close(60)

      write(*,*) "qfilter terminated normally."

      end program qfilter

      function gweight(x,dx,gy,gnum)
      implicit none
      integer maxsite
      parameter (maxsite=10000)
      integer i,j,k,gnum
      integer ixa,ixb
      integer pbead
      DOUBLE PRECISION x,gweight,dx,gy(0:gnum)
      DOUBLE PRECISION xa,xb,ya,yb,gmax,absx
      DOUBLE PRECISION loc,poc,relx,seg,sumu,sumb
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION xk(0:1000),xksq(0:1000)
      DOUBLE PRECISION aex,bex,egyik,masik


      common pbead

      PI=4.D0*DATAN(1.D0)

      absx=abs(x)

      gmax = dx*gnum

      if (absx .eq. gmax) then

          gweight=gy(gnum)

          return

      endif

      if (absx .gt. gmax) then

c    extrapolation

c real extrapolation

c           loc = gy(gnum-1) - (absx-dx)/2/pbead

c           poc = gy(gnum) - gmax/2/pbead

c           bex = log(loc/poc) / dx**2

c           bex=max(0.0d0,bex)

c           gweight = poc * exp(-bex*(absx-gmax)**2)

c           gweight = gweight + absx / 2 / pbead

c end of real extrapolation

c simple estimation

           gweight = absx / 2 / pbead

c end of simple estimation

      else

          ixa = int(absx/dx)
          ixb = ixa + 1
          xa= ixa * dx
          xb= ixb * dx

          ya= gy(ixa)
          yb= gy(ixb)

c linear interpolation:

          gweight = ((absx-xa)*yb+(xb-absx)*ya)/(xb-xa)

c          write(*,*) absx,ixa,ixb,xa,xb,ya,yb,gweight

      endif

      return
      end
