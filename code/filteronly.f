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
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: cx(:),cy(:),cz(:)
      character*120 xyzfile,outxyzfile,stratoms
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
  10    format (/,' Name of the input xyz file: ',$)
      read (*,20)  xyzfile
   20    format (a120)

      write (*,*) xyzfile

c
c     get the number of frames
c

      OPEN (50, file = xyzfile)
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

      write (*,80)
   80 format (/,' Timestep between frames',
     &              ' in picoseconds:  ',$)
   90 format (f20.0)
      read (*,90) timestep
      write (*,*) timestep

      write (*,81)
   81 format (/,' Temperature for filtration',
     &             ' in Kelvin:  ',$)
   91    format (f20.0)
         read (*,91) temp
         write (*,*) temp
      write(*,*) 'Calculation of the kernel function:'

      deltat=2*PI/(avogadro*planck/boltzmann/temp/timestep*1E11)

c Reading the kernel function

      INQUIRE(FILE=kernel_file_name, EXIST=kernel_file_exist)

      if (kernel_file_exist) then

        OPEN (51, file = kernel_file_name)

         j = 0
         DO
           READ (51,*, END=14) kernel(j)
           j = j + 1
         END DO
   14 CLOSE (51)

      else

          write(*,*) "kernel.dat file is not found. Program stops."
          stop

      endif

      ndec=j-1

      ablak=2*ndec+1

c Check the integral of the kernel function

      sum=kernel(0)/2.0d0

      do j=1,ndec

          sum=sum+kernel(j)

      end do

      sum=sum*2.0d0

      write(*,*) "The integral of the kernel function is: ",sum
      write(*,*) "before it is normalized to 1."
c Normalize the kernel function to 1.0

      do j=0,ndec

          write(*,*) kernel(j)

          kernel(j)=kernel(j)/sum

      end do

      write(*,*) 'Starting filtration'
c     End of the calculation of the kernel function



c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock))
      allocate (nam(natoms))
      allocate (cx(natoms))
      allocate (cy(natoms))
      allocate (cz(natoms))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      open (unit=50,file=xyzfile,status='old',action='read')  

      outxyzfile='filt_'//xyzfile

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
   15    format (a300)

         iline=iline+2

         do j=1,natoms
            iline=iline+1
            read(50,*,iostat=ierror) nam(j),cx(j),cy(j),cz(j)  

            if ( ierror .ne. 0 ) then

                write (*,*) "ERROR: Problem reading file ",xyzfile
                write (*,*) "in line ",iline
                write (*,*)  "Program stops."
            stop

            endif

         end do  
 
         ipos=mod(i-1,ablak)+1
         midpos=mod(i-1-ndec,ablak)+1
         do j = 1, natoms
            xtot(j,ipos)=cx(j)
            ytot(j,ipos)=cy(j)
            ztot(j,ipos)=cz(j)
         end do

         if (i .ge. ablak) then

            do m=1, natoms
               cx(m)=0.0d0
               cy(m)=0.0d0
               cz(m)=0.0d0
            end do
               
            do k=1,ablak
               tav=min(abs(midpos-k),ablak+k-midpos,ablak+midpos-k)
               do m=1, natoms
                  cx(m)=cx(m)+xtot(m,k)*kernel(tav)
                  cy(m)=cy(m)+ytot(m,k)*kernel(tav)
                  cz(m)=cz(m)+ztot(m,k)*kernel(tav)
               end do
            end do

c
c     write output
c

      write(*,*) 'Writing frame', i-2*ndec
      
            write(60,'(a)') trim(stratoms)

            write(60,'(a)') trim(title(i-ndec))

            do k=1,natoms
               write(60,25) nam(k),cx(k),cy(k),cz(k)  
   25    format (a3,3f11.5)
            end do  

         end if
         
      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)


      close(50)
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
