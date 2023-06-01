program make_restart
implicit none
integer                         :: i,j,ios,n
character(len=128)              :: fname,buf,camadir,expdir,mapname
!-map variables
real                            :: gsize,west, north, east, south ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated
! real,allocatable                :: rivsto(:,:),fldsto(:,:) ! put to restart file
real,allocatable                :: rivsto(:,:),fldsto(:,:),damsto(:,:) ! Modified by Youjiang, put to restart file

real,allocatable                :: elevtn(:,:)

real,allocatable                :: rivlen(:,:),rivwth(:,:),rivsto_max(:,:),rivhgt(:,:)
real,allocatable                :: fldhgt(:,:,:)
integer,allocatable             :: oceanmask(:,:),fldstage(:,:)
real,allocatable                :: grid_area(:,:)
real,allocatable                :: nextdst(:,:)
integer,allocatable             :: nextX(:,:),nextY(:,:)
real,parameter                  :: g=9.80665,dt=86400.,man=0.03,man2=0.10,pdstmth=10000.
character(len=8)                :: yyyymmdd,onedaybef,onedayaft
real                            :: dhgtpre        !! private
character(len=3)                :: num_name, cal
character(len=10)               :: loop

real,allocatable                :: fldfrac(:,:)

character(len=1)                :: loopchar
integer                         :: ens_num,k
integer                         :: cor ! for parameter corruption
character(len=3)                :: opt ! for CaMa-Flood option

real,allocatable                :: xa(:,:),rivdph(:,:),flddph(:,:)

real                            :: hgt,pre,Across

! how the restart file is made
!   rivsto,fldsto => recalculation
!   damsto ==> copy from CaMa_out
!   rivout,fldout => average of ensembles (extracted from restart file) ## not used
!   rivdpt,fldsto => average of ensembles (extracted from restart file) ## not used

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
read(buf,"(A)") mapname

call getarg(7,buf)
read(buf,*) ens_num ! number of ensemble members

call getarg(8,buf)
read(buf,*) num_name

call getarg(9,buf)
read(buf,"(A)") expdir

! call getarg(10,buf) ! not needed
! read(buf,*) cor ! for identifying the corruption 0: none, 1: rivhgt, 2: rivwth 
!                 ! 3:rivman 4: fldhgt, 5: rivhgt, rivwth, rivman, and fldhgt

call getarg(10,buf)
read(buf,"(A)") opt
!==
fname=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
open(11,file=fname,form='formatted')
read(11,*) lonpx
read(11,*) latpx
read(11,*) nflp
read(11,*) gsize
read(11,*) west
read(11,*) east
read(11,*) south
read(11,*) north
close(11)
!-------

if(trim(adjustl(loop))=="open") loopchar="C"
if(trim(adjustl(loop))=="assim") loopchar="A"

! program for making restart file
fname=trim(adjustl(expdir))//"/logout/restartError_"//yyyymmdd//loopchar//num_name//".log"
open(82,file=fname,status='replace')
write(82,*) "make_restart.f90 Errors"

! read assimilated xa
allocate(xa(lonpx,latpx))
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
! read CaMa-Flood parametes
allocate(rivlen(lonpx,latpx),rivwth(lonpx,latpx),rivhgt(lonpx,latpx),fldhgt(lonpx,latpx,nflp))
allocate(elevtn(lonpx,latpx),nextX(lonpx,latpx),nextY(lonpx,latpx),nextdst(lonpx,latpx),grid_area(lonpx,latpx))

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is 1.39e4
else
    write(*,*) "no file rivlen"
    write(82,*) "no file rivlen at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth_gwdlr.bin"
! if (cor==2 .OR. cor==5) then
!     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth_corrupt.bin"
! end if
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
!      fname=trim(adjustl(camadir))//"map/"//trim(mapname)//"/rivhgt_"//num_name//loopchar//".bin"
!    else
!      fname=trim(adjustl(camadir))//"map/"//trim(mapname)//"/rivhgt.bin"
!    end if
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt.bin"
! if (cor==1 .OR. cor==5) then
!     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt_corrupt.bin"
! end if
! ! !not needed for virtual experiments
! ! if (trim(cal)=="yes") then
! !     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt_Xudong.bin"
! ! elseif (trim(cal)=="corrupt") then
! !     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt_corrupt.bin"
! ! else
! !     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt.bin"
! ! end if
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivhgt
    ! ocean is -9999
else
    write(*,*) "no file rivhgt"
    write(82,*) "no file rivhgt at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/fldhgt.bin"
! if (cor==4 .OR. cor==5) then
!     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/fldhgt_corrupt.bin"
! end if
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx*10,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) fldhgt
    ! ocean is -9999
else
    write(*,*) "no file fldhgt"
    write(82,*) "no file fldhgt at:",fname
end if
close(34)

!    fname=trim(adjustl(camadir))//"map/"//trim(mapname)//"/elevtn.bin"
!    write(numch,'(i3.3)') num
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/elevtn.bin"
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
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nextxy.bin"
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
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nxtdst.bin"
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
!fname=trim(adjustl(camadir))//"map/"//trim(mapname)//"/grdare.bin"
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/ctmare.bin" ! after CaMa-Flood v3.9
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
! allocate
allocate(rivsto_max(lonpx,latpx),oceanmask(lonpx,latpx),fldstage(lonpx,latpx))
! allocate(fldfrac(lonpx,latpx),rivdph(lonpx,latpx),rivsto(lonpx,latpx),flddph(lonpx,latpx),fldsto(lonpx,latpx))
! Modified by Youjiang Shen
allocate(fldfrac(lonpx,latpx),rivdph(lonpx,latpx),rivsto(lonpx,latpx),damsto(lonpx,latpx),flddph(lonpx,latpx),fldsto(lonpx,latpx))


! calc river storage max
rivsto_max = rivlen*rivwth*rivhgt

! make ocean mask
! 0:inland 1:mouth 2:ocean
do i=1,lonpx
    do j=1,latpx
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
do i=1,lonpx
    do j=1,latpx
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
do i=1,lonpx
  do j=1,latpx
    if(fldstage(i,j)>0)then
      flddph(i,j) = rivdph(i,j) - rivhgt(i,j)
      !if (flddph(i,j)<1e-2) flddph(i,j) =0
    end if
  end do
end do

! calculate flood inundation percentage
do i=1,lonpx
  do j=1,latpx
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
do i=1,lonpx
    do j=1,latpx
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
do i=1,lonpx
    do j=1,latpx
        if(oceanmask(i,j)==2)then
            rivsto(i,j)=1e20
            fldsto(i,j)=1e20
        end if
    end do
end do

! =================================================
! Needed for CaMa-Flood v4.1 without reservoir assimilation --> may change in future
if (opt == "all" .OR. opt == "dam") then
    ! need to open CaMa_out restrat file to copy damsto
    ! open restart file at CaMa_out 
    print*, "make_restart.f90: open damsto"
    fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//loopchar//num_name//"/restart"//onedayaft//".bin"
    print*, "CaMa_out", fname
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        print* , "read damsto", fname
        read(34,rec=3) damsto  ! added by Youjiang. 
    else
        print*, "no file", fname
    end if
    close(34)
end if

! =================================================
! save and store output & restart file
! make restart file
! rivsto,fldsto,damsto, ... , ---> for CaMa-Flood V4.1.0
fname=trim(adjustl(expdir))//"/CaMa_in/restart/"//trim(adjustl(loop))//"/restart"//onedayaft//loopchar//num_name//".bin"
open(35,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(35,rec=1) rivsto
    write(35,rec=2) fldsto
    if (opt == "all" .OR. opt == "dam") then
        write(35,rec=3) damsto ! for CaMa-Flood v4.1
    end if
end if
close(35)
write(82,*) "done restart file at:",fname
print*, "---> restart ", fname
close(82)
! deallocate
deallocate(xa,rivlen,rivwth,rivhgt,fldhgt)
deallocate(elevtn,nextX,nextY,nextdst,grid_area)
deallocate(rivdph,rivsto,flddph,fldsto,damsto)
end program make_restart