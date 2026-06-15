program mean
!******************************************
! Calculate mean of a given variable
! inspired by Hansaki et al,. (2018)
! created by Menaka@IIS
!******************************************
implicit none
character(len=128)              :: fname,buf,camadir,expdir,mapname,outdir,tag !DAdir,patchdir,hydrowebdir,varname
character(len=8)                :: yyyymmdd,nxtyyyymmdd
character(len=3)                :: numch
integer                         :: syear, eyear, smon, emon, eday, N, ios
integer                         :: num, year, mon, day
integer                         :: isleap,leap,monthday,ensnum
!-map variables
real                            :: gsize,west, north, east, south ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated
integer                         :: nXX, nYY
real,allocatable                :: var(:,:,:)
!==================================
if ( iargc().lt.5 ) then
  write(*,*) 'Usage: mean'
  write(*,*) '      syear eyear CaMadir mapname expdir numch outdir tag N'
  stop
end if

call getarg(1,buf)
read(buf,*) syear ! start year

call getarg(2,buf)
read(buf,*) eyear ! end year

call getarg(3,buf)
read(buf,"(A)") camadir ! CaMadir

call getarg(4,buf)
read(buf,"(A)") mapname ! name of the map

call getarg(5,buf)
read(buf,"(A)") expdir

call getarg(6,buf)
read(buf,*) ensnum

call getarg(7,buf)
read(buf,"(A)") outdir

call getarg(8,buf)
read(buf,"(A)") tag

call getarg(9,buf)
read(buf,*) N

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
nXX=lonpx
nYY=latpx
!==
write(numch,'(i3.3)') ensnum
! print*,numch
!==
allocate(var(nXX,nYY,N))
!================
num=1
do year=syear,eyear
    isleap=leap(year)
    do mon=1, 12
        eday=monthday(mon,isleap)
        do day=1, eday
            write(yyyymmdd,'(i4.4,i2.2,i2.2)') year,mon,day
            fname=trim(expdir)//"/assim_out/ens_xa/open/"//trim(yyyymmdd)//"_"//trim(numch)//"_xa.bin"
            open(34,file=fname,form="unformatted",access="direct",recl=4*nXX*nYY,status="old",iostat=ios)
            if(ios==0)then
                read(34,rec=1) var(:,:,num)
            else
                write(*,*) "no: ",trim(fname)
            end if
            close(34)
            num=num+1
        end do
    end do
end do
!----------
! save file
!----------
!=======
! mean
fname=trim(adjustl(outdir))//"/mean_sfcelv_"//trim(tag)//"_"//trim(numch)//".bin"
write(*,*)trim(fname)
open(35,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(35,rec=1) sum(var,DIM=3)/(real(N)+1e-20)
else
    write(*,*) "not created", trim(fname)
end if
close(35)
!======
! standard deviation
fname=trim(adjustl(outdir))//"/std_sfcelv_"//trim(tag)//"_"//trim(numch)//".bin"
write(*,*)fname
open(35,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(35,rec=1) sqrt(sum(var**2,DIM=3)/(real(N)+1e-20) - (sum(var,DIM=3)/(real(N)+1e-20))**2)
else
    write(*,*) "not created", fname
end if
close(35)
! print*,numch,ensnum
deallocate(var)
end program mean
!*********************************
function leap(year)
  implicit none
  integer                        :: year,leap
  !real                           :: mod
  !--
  ! days of the year
  ! normal : 365
  ! leap   : 366
  leap=0
  if (mod(dble(year),4.0)   == 0) leap=1
  if (mod(dble(year),100.0) == 0) leap=0
  if (mod(dble(year),400.0) == 0) leap=1
  return
end function leap
!*********************************
function monthday(imon,leap)
  implicit none
  integer                        :: imon, leap, monthday
  ! ================================================
  ! Calculation for months except February
  ! ================================================
  if(imon.eq.4.or.imon.eq.6.or.imon.eq.9.or.imon.eq.11)then
      monthday=30
  else
      monthday=31
  end if
  ! ================================================
  ! Calculation for February
  ! ================================================
  if ( leap.eq.1)then
      if (imon.eq.2)then
          monthday=29
      end if
  else
      if (imon.eq.2)then
          monthday=28
      end if
  end if
  return
end function monthday
!*********************************