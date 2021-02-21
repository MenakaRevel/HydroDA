program mean
!******************************************
! Calculate mean of a given variable
! inspired by Hansaki et al,. (2018)
! created by Menaka@IIS
!******************************************
implicit none
character(len=128)              :: fname,buf,camadir,expdir,DAdir,patchdir,hydrowebdir,mapname,varname
character(len=8)                :: yyyymmdd,nxtyyyymmdd
!-map variables
real                            :: gsize,west, north, east, south ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated
!==================================
if(iargc().lt.5) then
  write(*,*) 'Usage: mean'
  write(*,*) '      yearmin yearmax yearout CaMadir mapname'
  stop
end if

call getarg(1,buf)
read(buf,*) yearmin ! minimum year

call getarg(2,buf)
read(buf,*) yearmax ! maximum year

call getarg(3,buf)
read(buf,*) yearout ! out year

call getarg(4,buf)
read(buf,"(A)") camadir ! CaMadir

call getarg(5,buf)
read(buf,*) mapname ! name of the map

call getarg(6,buf)
read(buf,"(A)") expdir

call getarg(7,buf)
read(buf,"(A)") DAdir

call getarg(8,buf)
read(buf,*) varname

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

end program mean
