program make_restart
implicit none
integer                         :: i,j,ios,n
integer,parameter               :: latpx=720,lonpx=1440
real,dimension(lonpx,latpx)     :: rivsto,fldsto ! put to restart file
real,dimension(lonpx,latpx)     :: elevtn
character*128                   :: fname,buf,camadir,expdir
real,dimension(lonpx,latpx)     :: rivlen,rivwth,rivsto_max,rivhgt
real,dimension(lonpx,latpx,10)  :: fldhgt
integer,dimension(lonpx,latpx)  :: oceanmask,fldstage
real,dimension(lonpx,latpx)     :: grid_area
real,dimension(lonpx,latpx)     :: nextdst
integer*4,dimension(lonpx,latpx):: nextX,nextY
real,parameter                  :: g=9.80665,dt=86400.,man=0.03,man2=0.10,pdstmth=10000.
character*8                     :: yyyymmdd,onedaybef,onedayaft
real                            :: dhgtpre        !! private
character*3                     :: num_name
character*10                    :: loop

real,dimension(lonpx,latpx)     :: fldfrac

character*1                     :: loopchar
integer                         :: ens_num,k

real,dimension(lonpx,latpx)     :: xa,rivdph,flddph

real                            :: hgt,pre,Across

! how the restart file is made
!   rivsto,fldsto => recalculation
!   rivout,fldout => average of ensembles (extracted from restart file)
!   rivdpt,fldsto => average of ensembles (extracted from restart file)

call getarg(1,buf)
read(buf,*) yyyymmdd

call getarg(2,buf)
read(buf,*) onedaybef

call getarg(3,buf)
read(buf,*) onedayaft

call getarg(4,buf)
read(buf,*) loop

call getarg(5,buf)
read(buf,"(A)") camadir

call getarg(6,buf)
read(buf,*) ens_num ! number of ensemble members

call getarg(7,buf)
read(buf,*) num_name

call getarg(8,buf)
read(buf,"(A)") expdir


if(trim(adjustl(loop))=="open") loopchar="C"
if(trim(adjustl(loop))=="assim") loopchar="A"

! program for making restart file
fname=trim(adjustl(expdir))//"/logout/restartError_"//yyyymmdd//loopchar//num_name//".log"
open(82,file=fname,status='replace')
write(82,*) "make_restart.f90 Errors"

! read assimilated xa
fname=trim(adjustl(expdir))//"/assim_out/ens_xa/"//trim(adjustl(loop))//"/"//yyyymmdd//"_"//num_name//"_xa.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) xa
else
    write(*,*) "no file at xa"
    write(82,*) "no file xa",fname
end if
close(34)

! read many parameters
fname=trim(adjustl(camadir))//"map/glb_15min/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is 1.39e4
else
    write(*,*) "no file rivlen"
    write(82,*) "no file rivlen at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"map/glb_15min/rivwth_gwdlr.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
    write(82,*) "no file rivwth at:",fname 
end if
close(34)

!    if ((loopchar == "C") .or. (loopchar == "A")) then 
!      fname=trim(adjustl(camadir))//"map/glb_15min/rivhgt_"//num_name//loopchar//".bin"
!    else
!      fname=trim(adjustl(camadir))//"map/glb_15min/rivhgt.bin"
!    end if
fname=trim(adjustl(camadir))//"map/glb_15min/rivhgt.bin" 
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivhgt
    ! ocean is -9999
else
    write(*,*) "no file rivhgt"
    write(82,*) "no file rivhgt at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"map/glb_15min/fldhgt.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx*10,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) fldhgt
    ! ocean is -9999
else
    write(*,*) "no file fldhgt"
    write(82,*) "no file fldhgt at:",fname
end if
close(34)

!    fname=trim(adjustl(camadir))//"map/glb_15min/elevtn.bin"
!    write(numch,'(i3.3)') num
fname=trim(adjustl(camadir))//"map/glb_15min/elevtn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) elevtn
    ! ocean is -9999
else
    write(*,*) "no file elevtn"
    write(82,*) "no file elevtn at:",fname
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"map/glb_15min/nextxy.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextX
    read(34,rec=2) nextY
else
    write(*,*) "no file nextXY"
    write(82,*) "no file nextXY at:",fname
end if
close(34)

! read distance to next grid
fname=trim(adjustl(camadir))//"map/glb_15min/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
    ! 海は
else
    write(*,*) "no file nextdst"
    write(82,*) "no file nextdst at:",fname
end if
close(34)


! read grid area
!fname=trim(adjustl(camadir))//"map/glb_15min/grdare.bin"
fname=trim(adjustl(camadir))//"map/glb_15min/ctmare.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) grid_area
    ! 海は
else
    write(*,*) "no file grid_area"
    write(82,*) "no file grarea at:",fname
end if
close(34)

! =======================================================
! calc river storage max
rivsto_max = rivlen*rivwth*rivhgt

! make ocean mask
! 0:inland 1:mouth 2:ocean
do i=1,1440
    do j=1,720
        if(nextX(i,j)==-9999)then
            oceanmask(i,j)=2
        else if(nextX(i,j)==-9)then ! ocean
            oceanmask(i,j)=1
        else if(nextX(i,j)==-10)then ! inland termination
            oceanmask(i,j)=1
        else
            oceanmask(i,j)=0
        end if
    end do
end do

! =======================================================
! NEW judgement of flood stage from WSE
! 0: xa <= max  1: xa > max
! WSE = elevtn - rivhgt + rivdph
! hence, max = elevtn
fldstage=0
fldfrac=0
do i=1,1440
    do j=1,720
        if(oceanmask(i,j)<2)then
            if(xa(i,j)>elevtn(i,j))then
                fldstage(i,j)=1
            end if
        end if
    end do
end do


! calculate rivdph
rivdph = xa - elevtn + rivhgt

! calculate river storage
rivsto = rivdph * rivlen * rivwth

! calculate flddph
flddph = 0
do i=1,1440
  do j=1,720
    if(fldstage(i,j)>0)then
      flddph(i,j) = rivdph(i,j) - rivhgt(i,j)
      !if (flddph(i,j)<1e-2) flddph(i,j) =0
    end if
  end do
end do

! calculate flood inundation percentage
do i=1,1440
  do j=1,720
    if(fldstage(i,j)==1)then
      do n=1,10
        if(fldhgt(i,j,n)>flddph(i,j)) exit
        fldstage(i,j) = fldstage(i,j) + 1
      end do
    end if
    ! fldstage(i,j) == 0: no flood
    ! fldstage(i,j) 1~10: flood stage |  1  |  2  |  ...  |  10  |
    !                          fldhgt ↑0m  ↑[1] ↑[2] ..↑[9]  ↑[10]
    ! fldstage(i,j) ==11: over fldmax
  end do
end do

! calculate fldsto
fldsto = 0
do i=1,1440
  do j=1,720
    ! calculate flood storage until fldstage(i,j)-1
    if(fldstage(i,j)>1)then
      dhgtpre = 0
      do n=1,fldstage(i,j)-1
        fldsto(i,j) = fldsto(i,j) + grid_area(i,j)*0.1*(real(n)-0.5)*(fldhgt(i,j,n)-dhgtpre)
        dhgtpre = fldhgt(i,j,n)
        fldfrac(i,j) = fldfrac(i,j) + 0.1
      end do
    end if

    ! calculate flood storage at fldstage(i,j)
    if(fldstage(i,j)>0 .and. fldstage(i,j)<11)then
      if(fldstage(i,j)==1)then
        dhgtpre = 0
      end if
      if(fldstage(i,j)>1)then
        dhgtpre = fldhgt(i,j,fldstage(i,j)-1)
      end if
      n = fldstage(i,j)
      k = fldstage(i,j)-1
      hgt = fldhgt(i,j,n)
      pre = dhgtpre
 
      !Across = grid_area(i,j)*(k*(hgt-flddph(i,j))+(K+1)*(flddph(i,j)-pre))/(hgt-pre+1e-20)*0.1
      Across = 0.1*grid_area(i,j)*(k + ((flddph(i,j)-pre)/(hgt-pre+1e-20))*1.0)    !(k*(hgt-flddph(i,j))+(k+1)*(flddph(i,j)-pre))/(hgt-pre+1e-20)
      fldfrac(i,j) = fldfrac(i,j) + 0.1/(hgt-pre+1e-20)*(0*(hgt-flddph(i,j))+1*(flddph(i,j)-pre))

      fldsto(i,j) = fldsto(i,j) + ((Across+0.1*grid_area(i,j)*k)*0.5)*(flddph(i,j)-pre)

    end if
    if (fldstage(i,j)==11)then
      pre = fldhgt(i,j,10)
      fldsto(i,j) = fldsto(i,j) + grid_area(i,j)*(flddph(i,j)-pre)
      fldfrac(i,j) = 1.0
    end if

  end do
end do

! ===================================================
! arrangement for ocean mask
do i=1,1440
    do j=1,720
        if(oceanmask(i,j)==2)then
            rivsto(i,j)=1e20
            fldsto(i,j)=1e20
        end if
    end do
end do

! =================================================
! save and store output & restart file

! make restart file
! rivsto,fldsto,rivout,fldout,rivdpt_pre,fldsto_pre
fname=trim(adjustl(expdir))//"/CaMa_in/restart/"//trim(adjustl(loop))//"/restart"//onedayaft//loopchar//num_name//".bin"
open(35,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(35,rec=1) rivsto
    write(35,rec=2) fldsto
end if
close(35)
write(82,*) "done restart file at:",fname
close(82)
end program make_restart
