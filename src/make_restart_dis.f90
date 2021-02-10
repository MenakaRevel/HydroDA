program make_restart_dis
implicit none
integer                         :: i,j,ios,n,ix,iy,iix,iiy
character(len=128)              :: fname,buf,camadir,expdir,mapname
!-map variables
real                            :: gsize,west, north, east, south ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated

real,allocatable                :: rivsto(:,:),fldsto(:,:), sto(:,:) ! put to restart file
real,allocatable                :: prerivsto(:,:),prefldsto(:,:), presto(:,:) ! read restart file
real,allocatable                :: elevtn(:,:)

real,allocatable                :: rivlen(:,:),rivwth(:,:),rivsto_max(:,:),rivhgt(:,:)
real,allocatable                :: fldhgt(:,:,:)
integer,allocatable             :: oceanmask(:,:),fldstage(:,:)
real,allocatable                :: grid_area(:,:)
real,allocatable                :: nextdst(:,:)
integer,allocatable             :: nextX(:,:),nextY(:,:),rivseq(:,:)
real,parameter                  :: g=9.80665,dt=86400.,man=0.03,man2=0.10,pdstmth=10000.
character(len=8)                :: yyyymmdd,onedaybef,onedayaft
real                            :: dhgtpre        !! private
character(len=3)                :: num_name
character(len=10)               :: loop

real,allocatable                :: fldfrac(:,:)

character(len=1)                :: loopchar
integer                         :: ens_num,k

real,allocatable                :: xa(:,:),xf(:,:),rivdph(:,:),flddph(:,:)
integer,allocatable             :: xlist(:),ylist(:)
real                            :: hgt,pre,Across
real                            :: disin,disout,predisin,predisout

real                            :: stopre, wthpre, dphpre, wthinc
real                            :: stonow, wthnow
real,allocatable                :: fldstomax(:), fldgrd(:)

! how the restart file is made
!   rivsto,fldsto => recalculation

! Assume outflow is the is the instantaneous flow at the final time of the
! integration window. Assume stroage differnce due dischare differnce in 
! in flow discharge to the gird (ix,iy) and out flow discharge from (ix,iy).

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
! this should be discharge
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

! read discahrge before assimilation
allocate(xf(lonpx,latpx))
fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//trim(loopchar)//num_name//"/outflw"//yyyymmdd(1:4)//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) xf
else
    write(*,*) "no file at xf"
    write(82,*) "no file xf",fname
end if
close(34)

! read prestorage 
allocate(prerivsto(lonpx,latpx),prefldsto(lonpx,latpx))
fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//trim(loopchar)//num_name//"/restart"//yyyymmdd//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) prerivsto
    read(34,rec=2) prefldsto
else
    write(*,*) "no file at xf"
    write(82,*) "no file xf",fname
end if
close(34)


! read many parameters
! read CaMa-Flood parametes
allocate(rivlen(lonpx,latpx),rivwth(lonpx,latpx),rivhgt(lonpx,latpx),fldhgt(lonpx,latpx,nflp),rivseq(lonpx,latpx))
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
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/ctmare.bin" ! after CaMa v3.9
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) grid_area
    ! 海は
else
    write(*,*) "no file grid_area"
    write(82,*) "no file grarea at:",fname
end if
close(34)

! read river sequence
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivseq.bin" 
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivseq
    ! 海は
else
    write(*,*) "no file rivseq"
    write(82,*) "no file rivseq at:",fname
end if
close(34)

! =======================================================
! allocate
allocate(presto(lonpx,latpx),sto(lonpx,latpx),rivsto_max(lonpx,latpx),oceanmask(lonpx,latpx),fldstage(lonpx,latpx))
allocate(fldfrac(lonpx,latpx),rivdph(lonpx,latpx),rivsto(lonpx,latpx),flddph(lonpx,latpx),fldsto(lonpx,latpx))
allocate(xlist(8),ylist(8),fldstomax(nflp),fldgrd(nflp))
presto= rivsto+ fldsto
sto= presto
! make calculation
do ix=1,lonpx
    do iy=1,latpx
        if (nextX(i,j) == -9999.0) then
            cycle
        end if
        !----
        ! upstram inflow
        xlist=-9999
        ylist=-9999
        predisin=-9999.0
        predisout=-9999.0
        disin=-9999.0
        disout=-9999.0
        if (rivseq(ix,iy)==1) then
            print*, "rivseq : 1"
            predisin=0.0*dt
            disin=0.0*dt
        elseif (rivseq(ix,iy)>1) then
            call upgrids(ix,iy,lonpx,latpx,nextX,nextY,xlist,ylist,k)
            predisin=0.0
            do i=1,k
                iix=xlist(i)
                iiy=ylist(i)
                predisin=predisin+xf(iix,iiy)*dt
                disin=disin+xa(iix,iiy)*dt
            end do
        end if
        predisout=xf(ix,iy)*dt
        disout=xa(ix,iy)*dt
        sto(ix,iy) = presto(ix,iy) + disin - predisin - (disout - predisout)
    end do
end do        

! calc river storage max
rivsto_max = rivlen*rivwth*rivhgt
do ix=1,lonpx
    do iy=1,latpx
        if (nextX(i,j) == -9999.0) then
            cycle
        end if
        if (sto(ix,iy) > rivsto_max(ix,iy)) then
            i=1
            stopre = rivsto_max(ix,iy)
            wthpre = rivwth(ix,iy)
            dphpre = 0.0
            wthinc = (grid_area(ix,iy)/rivlen(ix,iy))*(1.0/real(nflp))
            call calc_fldstg(rivsto_max(ix,iy),grid_area(ix,iy),rivlen(ix,iy),rivwth(ix,iy),nflp,fldhgt(ix,iy,:),fldstomax,fldgrd)
            do while (sto(ix,iy) > fldstomax(i) .and. i<=nflp)
                stopre = fldstomax(i)
                wthpre = wthpre + wthinc
                dphpre = dphpre + fldhgt(ix,iy,i)
                i=i+1
                if (i>nflp) exit
            end do
            if (i>nflp) then
                stonow = sto(ix,iy) - stopre
                wthnow = 0.0
                flddph(ix,iy) = dphpre + (stonow/(wthpre * rivlen(ix,iy)))
            else
                stonow = sto(ix,iy) - stopre
                wthnow = -wthpre + ( wthpre**2. + 2.d0 * stonow * rivlen(ix,iy)**(-1.) * fldgrd(i)**(-1.) )**0.5
                flddph(ix,iy) = dphpre + fldgrd(i) * wthnow
            end if
            rivsto(ix,iy)=rivsto_max(ix,iy) + rivlen (ix,iy) * rivwth(ix,iy) * flddph(ix,iy)
            fldsto(ix,iy)=max(sto(ix,iy)-rivsto(ix,iy),0.0)
        else
            flddph(ix,iy)=0.0
            rivsto(ix,iy)=sto(ix,iy)
            fldsto(ix,iy)=0.0
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
!deallocate
deallocate(xa,rivlen,rivwth,rivhgt,fldhgt)
deallocate(elevtn,nextX,nextY,nextdst,grid_area)
deallocate(rivdph,rivsto,flddph,fldsto)
end program make_restart_dis  
!=======================================================
subroutine upgrids(i,j,nx,ny,nextX,nextY,xlist,ylist,k)
implicit none
! calculate nearst upstream grids
integer                           :: i,j,nx,ny
integer,dimension(nx,ny)          :: nextX,nextY
integer,dimension(8)              :: xlist,ylist
integer                           :: k

integer                           :: ix,iy,dd

dd=1
do ix=i-dd,i+dd
    do iy=j-dd,j+dd
        if (nextX(ix,iy) == i .and. nextY(ix,iy)==j) then
            xlist(k)=ix
            ylist(k)=iy 
            k=k+1
        end if
    end do
end do
k=k-1
return
end subroutine upgrids
!=======================================================
subroutine calc_fldstg(rivstomax,grarea,rivlen,rivwth,nflp,fldhgt,fldstomax,fldgrd)
implicit none
!* local variables
integer              :: nflp
real                 :: rivstomax,grarea,rivlen,rivwth
real,dimension(nflp) :: fldhgt
real,dimension(nflp) :: fldstomax, fldgrd

integer              :: i
real                 :: dstonow, dstopre, dhgtpre, dwthinc
!================================================
fldstomax(:) = 0.d0
fldgrd(:) = 0.d0
dstopre = rivstomax
dhgtpre = 0.d0
dwthinc = grarea * rivlen**(-1.) * (1.0/real(nflp))
do i=1, nflp
    dstonow      = rivlen * ( rivwth + dwthinc*(dble(i)-0.5) ) * (fldhgt(i)-dhgtpre)
    fldstomax(i) = dstopre + dstonow
    fldgrd(i)    = (fldhgt(i)-dhgtpre) * dwthinc**(-1.)
    dstopre      = fldstomax(i)
    dhgtpre      = fldhgt(i)
end do
return
!
end subroutine calc_fldstg