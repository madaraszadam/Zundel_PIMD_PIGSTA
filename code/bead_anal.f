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
      integer i,j,k,l,m,ipos,tav,midpos,nsp,glim,nstmax
      integer ierror,iline
      integer pbead            ! P: number of replacs or beads
      integer nspacx,nspacy,nspacz
      DOUBLE PRECISION rscal,mass(3),totmass,xcom,ycom,zcom
      DOUBLE PRECISION pbcx,pbcy,pbcz
      DOUBLE PRECISION hpbcx,hpbcy,hpbcz
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: cx(:,:),cy(:,:),cz(:,:)
      DOUBLE PRECISION, allocatable :: dd(:,:,:)
      character*120 basename,xyzfile,outxyzfile,stratoms,stri
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

      DOUBLE PRECISION roh1,roh2,roo,rgyr
      DOUBLE PRECISION sqoh1,sqoh2,sqoo
      DOUBLE PRECISION alfa,spt,gyrO,gyrHs,gyrHd

      LOGICAL :: gx_file_exist,kernel_file_exist

      INTEGER, allocatable :: atyp(:)

      character*120 gx_file_name,kernel_file_name

      INTEGER outgnum
      INTEGER fnum

      fnum=500

      rscal=1.0d0

      mass(1)=15.999
      mass(2)=1.0080
      mass(3)=1.0080

      totmass=mass(1)+mass(2)+mass(3)


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

      kernel(0)=1.0d0

      ndec=0

      ablak=2*ndec+1

c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock))
      allocate (nam(natoms))
      allocate (atyp(natoms))
      allocate (cx(natoms,pbead))
      allocate (cy(natoms,pbead))
      allocate (cz(natoms,pbead))

      allocate (dd(natoms,natoms,pbead))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      do i=1,pbead

        write (stri,'(I0)') i 
        xyzfile = trim(basename)//trim(stri)//"-1.xyz"

        write(*,*) "reading file ", xyzfile

        open (fnum+i,file=xyzfile,status='old',action='read')  

      end do

      open (unit=60,file="rOHd.txt",status='unknown',action='write')  
      open (unit=61,file="alfa.txt",status='unknown',action='write')
      open (unit=62,file="delta.txt",status='unknown',action='write')
      open (unit=63,file="gyrO.txt",status='unknown',action='write')
      open (unit=64,file="gyrHd.txt",status='unknown',action='write')
      open (unit=65,file="gyrHs.txt",status='unknown',action='write')


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
           read(fnum+k,15) title(i)
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



         do j=1,natoms

c determine atom types
c atyp means atom type
c atyp=1 means oxygen
c atyp=2 means dangling hydrogen 
c atyp=3 means shared proton

           if(trim(nam(j)).eq."O") then

c             write(*,*) "O atom"

             atyp(j)=1

           else

             if(trim(nam(j)).eq."H") then

               atyp(j)=2

             else

               atyp(j)=0

             endif

           endif

         end do

c calculate the distance matrix for O, H atom pairs
 
         do j=1,natoms

           do k=1,natoms

             if((k.ne.j).and.(atyp(j).eq.1)) then
c             if((atyp(j).eq.1).and.(atyp(k).eq.2)) then

               do m=1,pbead

                 dd(j,k,m) = (cx(j,m)-cx(k,m))**2

                 dd(j,k,m) =  dd(j,k,m) + (cy(j,m)-cy(k,m))**2

                 dd(j,k,m) =  dd(j,k,m) + (cz(j,m)-cz(k,m))**2

                 dd(k,j,m) = dd(j,k,m)

c                 write(*,*) dsqrt(dd(j,k,m))

               end do

             endif

           end do

         end do

c let us find the shared proton(s)

         oseg=maxval(dd)*2.0d0*pbead

         do j=1,natoms

           do k=j+1,natoms

             if((atyp(j).eq.1).and.(atyp(k).eq.1)) then

               do l=1,natoms

                if(atyp(l).ne.1) then

                 seg=0.0d0

                 do m=1,pbead

                   seg=seg+dd(j,l,m)+dd(k,l,m)

                 end do

                 if (seg.lt.oseg) then

                   nsp=l

                   oseg=seg

                 endif

                endif

               end do

c               write(*,*) "index of the shared proton: ", nsp

               atyp(nsp)=3

               do m=1,pbead

                sqoh1=dd(j,nsp,m)
                sqoh2=dd(k,nsp,m)
                sqoo=dd(j,k,m)

                roh1=dsqrt(sqoh1)
                roh2=dsqrt(sqoh2)
                roo=dsqrt(sqoo)

                spt=roh1-roh2

                write(62,*) spt

c                write(*,*) "shared proton coordinate: ", spt

c                write(*,*) roh1, roh2, roo

                alfa=acos((sqoh1+sqoh2-sqoo)/roh1/roh2/2.0d0)

                alfa=alfa/PI*180.0d0

                write(61,*) alfa

               end do

             endif

           end do

         end do

c calculate the radius of gyration

         do j=1,natoms

           rgyr=0.0d0

           xcom=0.0d0
           ycom=0.0d0
           zcom=0.0d0

           oseg=maxval(dd)*pbead

           do m=1,pbead

             xcom = xcom + cx(j,m)
             ycom = ycom + cy(j,m)
             zcom = zcom + cz(j,m)

           end do

           xcom = xcom / dble(pbead)
           ycom = ycom / dble(pbead)
           zcom = zcom / dble(pbead)

           do m=1,pbead

             rgyr = rgyr + (xcom - cx(j,m))**2
             rgyr = rgyr + (ycom - cy(j,m))**2
             rgyr = rgyr + (zcom - cz(j,m))**2

           end do

           rgyr = dsqrt(rgyr/dble(pbead))

           if(atyp(j).eq.1) then

             write(63,*) rgyr

           else

             if(atyp(j).eq.2) then

               write(64,*) rgyr

               do k=1,natoms

                 if(atyp(k).eq.1) then

                   seg=0.0d0

                   do m=1,pbead

                     seg=seg+dd(j,k,m)

                   end do

                   if (seg.lt.oseg) then

                     nsp=k

                     oseg=seg

                   endif


                 endif

               end do             

               do m=1,pbead

                 write(60,*) dsqrt(dd(j,nsp,m))

               end do

             else

               write(65,*) rgyr

             endif

           endif

         end do


c
c     write output
c

      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (atyp)
      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (dd)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)

      close(60)
      close(61)
      close(62)
      close(63)
      close(64)
      close(65)

      close(fnum)

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
