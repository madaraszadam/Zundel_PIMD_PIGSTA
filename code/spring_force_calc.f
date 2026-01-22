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

      DOUBLE PRECISION boltzmann,hbar,kjmol,abohr

      parameter (boltzmann=1.3806504d-23)        ! [J/K]
      parameter (hbar=1.05457162825177E-34)      ! [J*s]
      parameter (kjmol=2.62549961709828d+03)     ! a.u. to kJ/mol conversion factor
      parameter (abohr=5.29177208590000d-01)     ! Angström to Bohr

      integer gnum        ! number of table for the weight function
      integer ndec,ablak,nblock,natoms
      integer i,j,k,l,m,ipos,tav,midpos,nsp,glim,nstmax
      integer ierror,iline
      integer pbead            ! P: number of replacs or beads
      integer nspacx,nspacy,nspacz
      DOUBLE PRECISION rscal,totmass,xcom,ycom,zcom
      DOUBLE PRECISION pbcx,pbcy,pbcz
      DOUBLE PRECISION hpbcx,hpbcy,hpbcz
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: spm(:),cx(:,:),cy(:,:),cz(:,:)
      DOUBLE PRECISION, allocatable :: dd(:,:,:)
      character*120 basename,xyzfile,outxyzfile,stratoms,stri
      character(300), allocatable :: title(:,:)  
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

      DOUBLE PRECISION roh1,roh2,roo,rgyr
      DOUBLE PRECISION sqoh1,sqoh2,sqoo
      DOUBLE PRECISION alfa,spt,gyrO,gyrHs,gyrHd

      DOUBLE PRECISION spk ! spring constant

      DOUBLE PRECISION fx,fy,fz

      INTEGER nbp,nbn ! index of the previous bead and the next

      LOGICAL :: gx_file_exist,kernel_file_exist

      INTEGER, allocatable :: atyp(:)

      character*120 gx_file_name,kernel_file_name

      INTEGER outgnum
      INTEGER fnum,spnum

      fnum=500

      rscal=1.0d0

      kernel_file_name="kernel.dat"

      PI=4.D0*DATAN(1.D0)

c
c     get the base name of user specified input structures
c
      write (*,10)
  10    format (/,' Name of the input xyz file: ',$)
      read (*,20)  basename
   20    format (a120)

      write (*,*) basename



c
c     get the number of frames
c

      xyzfile=trim(basename)//"1-1.xyz"

      OPEN (fnum, file = xyzfile)
      read(fnum,'(i10)',iostat=ierror) natoms

      if ( ierror .ne. 0 ) then

            write (*,*) "ERROR: Incorrect xyz file format. ",
     +          "The first line is not an integer ",
     +          "that should mean the number of atoms. Program stops."
            stop

      endif

      i = 1
      DO
       READ (fnum,*, END=11)
       i = i + 1
      END DO
   11 CLOSE (fnum)

      if (mod(i,(natoms+2)) .ne. 0) then

            write (*,*) "WARNING: The number of lines does not ",
     +          "match the number of atoms."
            stop

      endif 

      nblock=i/(natoms+2)

      write(*,*) 'Number of atoms: ',natoms
      write(*,*) 'Number of input frames: ',nblock

      write (*,80)
   80 format (/,' Number of beads: ',$)
   90 format (f20.0)
      read (*,'(i10)',iostat=ierror) pbead 
      write (*,*) pbead

      write (*,81)
   81 format (/,' Temperature for filtration',
     &             ' in Kelvin:  ',$)
   91    format (f20.0)
         read (*,91) temp
         write (*,*) temp

      kernel(0)=1.0d0

      ndec=0

      ablak=2*ndec+1

c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock,pbead))
      allocate (nam(natoms))
      allocate (atyp(natoms))
      allocate (spm(natoms))
      allocate (cx(natoms,pbead))
      allocate (cy(natoms,pbead))
      allocate (cz(natoms,pbead))

      allocate (dd(natoms,natoms,pbead))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      spnum=fnum+pbead+1


c calculate the spring constant

      spk=(pbead*temp*boltzmann/hbar)**2/1.0d23/kjmol*abohr

c      spk=spk*1.0079d-3

c      spk=(pbead*temp*boltzmann/hbar)**2/1.0d23/kjmol

c      spk=spk*1.0079d-3

c      write(*,*) "spring constant: ", spk," a.u./Angstrom/Bohr"

c      stop

      do i=1,pbead

        write (stri,'(I0)') i 

        xyzfile = trim(basename)//trim(stri)//"-1.xyz"

        outxyzfile = "spf_"//xyzfile

        write(*,*) "reading file ", xyzfile

        open (fnum+i,file=xyzfile,status='old',action='read')  

        open (spnum+i,file=outxyzfile,status='unknown',action='write')

      end do

      write(*,*) 'Name of the outputfile: ', outxyzfile
      write(*,*) 'Number of filtered frames: ', nblock-ablak+1

      iline=0

c
c     cycle over all pairs of snapshot frame blocks
c
      do i = 1, nblock

c      do i = 1, 1

         do k=1,pbead

           read(fnum+k,'(a)') stratoms
           read(fnum+k,15) title(i,k)
   15    format (a300)

           iline=iline+2

           do j=1,natoms
              iline=iline+1
              read(fnum+k,*,iostat=ierror) nam(j),
     +                                     cx(j,k),cy(j,k),cz(j,k)

              if ( ierror .ne. 0 ) then

                  write (*,*) "ERROR: Problem reading file ",xyzfile
                  write (*,*) "in line ",iline
                  write (*,*)  "Program stops."
              stop

              endif

           end do  

         end do

       if (i.eq.1) then

c determine the masses of the atoms

         do j=1,natoms

           if(trim(nam(j)).eq."O") then

             spm(j)=spk*15.9994d-3

           else

             if(trim(nam(j)).eq."H") then

               spm(j)=spk*1.0079d-3

             else

               if(trim(nam(j)).eq."C") then

                 spm(j)=spk*12.0d-3

               else

              write(*,*) "Atom ",nam(j)," without mass, program stops."

                 stop

               endif

             endif

           endif

         end do

       endif

c
c     write output
c

        do k=1,pbead

          write(*,*) 'Writing frame', i
      
          write(spnum+k,'(a)') trim(stratoms)

          write(spnum+k,'(a)') trim(title(i,k))


          if (k.eq.1) then

            nbp=pbead

          else

            nbp=k-1

          endif

          if (k.eq.pbead) then

            nbn=1

          else

            nbn=k+1

          endif


          do j=1,natoms

            fx=spm(j)*(cx(j,nbn)+cx(j,nbp)-2*cx(j,k))
            fy=spm(j)*(cy(j,nbn)+cy(j,nbp)-2*cy(j,k))
            fz=spm(j)*(cz(j,nbn)+cz(j,nbp)-2*cz(j,k))

            write(spnum+k,25) nam(j),fx,fy,fz

   25    format (a3,3f16.10)
          end do  

        end do

      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (atyp)
      deallocate (spm)
      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (dd)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)

      do i=1,pbead

        close(fnum+i)

        close(spnum+i)

      end do

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
