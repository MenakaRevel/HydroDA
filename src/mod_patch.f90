module patch
!====================================================================================
! purpose: Procedures realted to local patches: reading local patch, etc
!
! 
!====================================================================================
! created by Menaka
! Menaka@IIS 2023
!====================================================================================
!$ use omp_lib
use common

implicit none

public

contains
!************************************************************************************
subroutine read_elp(fname,num,xlist,ylist,wgt)
!=======================================================================
! reading emperical local patch file
! input 
!  fname - file name
!  num   - nuumber of pixels in the local patch
! output
!  xlist - list of x corrdinates
!  ylist - list of y corrdinates
!  wgt   - list of observation localzation weights
! emperical local patch file includes:
! 1 - x location
! 2 - y location
! 3 - observvation localization weight (using gaussian weight)
!=======================================================================
implicit none
!--in
character(len=128),intent(in)      :: fname
integer,intent(in)                 :: num
!--out
integer,intent(out)                :: xlist(num),ylist(num)
real,intent(out)                   :: wgt(num)
!--
integer                            :: ios
! initialize
xlist=0
ylist=0
wgt=0
open(34,file=fname,status='old',access='sequential',form='formatted',action='read',iostat=ios)!
if(ios/=0)then
    print*, "no local patch file", fname
    ! write(82,*) "no local patch file at:", fname
    ! write(78,*) "no local patch file at:", fname
    goto 1090
end if
1000 continue
    read(34,*,end=1090) xlist,ylist,wgt
    goto 1000
1090 continue
close(34)
return
end subroutine read_elp
!************************************************************************************
subroutine read_wgt(fname,nx,ny,weightage)
!$ use omp_lib    
implicit none
!--in
character(len=128),intent(in)      :: fname
integer,intent(in)                 :: nx,ny
!--out
real,dimension(nx,ny),intent(out)  :: weightage
!--
integer                            :: ios
open(34,file=fname,form="unformatted",access="direct",recl=4*ny*nx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) weightage
else
    write(*,*) "no weightage", fname
end if
close(34)
return
end subroutine read_wgt
!************************************************************************************
function Gauss_wt(lag)
implicit none
real                                :: lag,Gauss_wt
real,parameter                      :: sigma=1000.0 !1000 km 
!---
Gauss_wt=exp(-(lag**2.0/(2.0*sigma**2.0)))  
!---
return
!---
end function Gauss_wt
!************************************************************************************
subroutine assign_local_patch(countnumber,targetpixel,patch_size,patch_start,patch_end,target_pixel,countnum)
!=======================================================================
! Assining local patches
implicit none
!--in
integer,intent(in)                 :: countnumber,targetpixel,patch_size
!--out
integer,intent(out)                :: patch_start,patch_end,target_pixel,countnum
!--
if (patch_size == 0) then ! for zero local patch ***Only target pixel is used
    patch_start=targetpixel
    patch_end=targetpixel
    target_pixel=1
    countnum=1
else
    patch_start=1
    patch_end=countnumber
    target_pixel=targetpixel
    countnum=countnumber
end if
return
end subroutine assign_local_patch
!************************************************************************************
end module patch