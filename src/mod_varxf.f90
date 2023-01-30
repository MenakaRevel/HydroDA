module varxf
!====================================================================================
! purpose: Procedures related to getting modeled variable
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
subroutine local_xf(globalx,xlist,ylist,countnum,patch_start,patch_end,nx,ny,ne,xf)
!=======================================================================
! assing local prognostic variable
!=======================================================================
implicit none
!--in
integer,intent(in)                             :: countnum,patch_start,patch_end,nx,ny,ne !,conflag
integer,intent(in)                             :: xlist(countnum),ylist(countnum)
real(r_size),intent(in)                        :: globalx(nx,ny,ne)!,meanglobalx(nx,ny,ne),stdglobalx(nx,ny,ne)
!--out
!integer,intent(out)                            :: local_sat(countnum)
real(r_size),intent(out)                       :: xf(countnum,ne) !,local_err(countnum)
!--
integer                                        :: i,j,i_m,j_m,num
xf=0
!print*,"L538: read model forcasts"
j=1
do i=patch_start,patch_end
    i_m=xlist(i)
    j_m=ylist(i)
    do num=1, ne
        xf(j,num)=globalx(i_m,j_m,num)
        ! if (conflag == 1) then
        !     xf(j,num)=globalx(i_m,j_m,num)
        ! else if (conflag == 2) then
        !     xf(j,num)=globalx(i_m,j_m,num)-meanglobalx(i_m,j_m,num) !meanglobaltrue(i_m,j_m)
        ! else if (conflag == 3) then
        !     ! xf(j,num)=(globalx(i_m,j_m,num)-meanglobaltrue(i_m,j_m))/(stdglobaltrue(i_m,j_m)+1.0e-20)
        !     xf(j,num)=(globalx(i_m,j_m,num)-meanglobalx(i_m,j_m,num))/(stdglobalx(i_m,j_m,num)+1.0e-20)
        ! else if (conflag == 4) then
        !     xf(j,num)=log10(globalx(i_m,j_m,num))
        ! end if
    end do
    j=j+1
end do
return
end subroutine local_xf
!************************************************************************************
subroutine get_ensemble_mean(x,countnum,ne,xm)
!=======================================================================
! get ensemble mean of x
!=======================================================================
!--in
integer,intent(in)                             :: ne,countnum
real(r_size),intent(in)                        :: x(countnum,ne)
!--out
real(r_size),intent(out)                       :: xm(countnum)
!--
integer                                        :: i
!---
do i=1,countnum
    xm(i)=sum(x(i,:))/(1e-20+real(ne))
end do
return
end subroutine get_ensemble_mean
!************************************************************************************
subroutine get_ensemble_diff(x,xm,countnum,ne,E)
!=======================================================================
! get ensemble mean differnce of x
!=======================================================================
!--in
integer,intent(in)                             :: ne,countnum
real(r_size),intent(in)                        :: x(countnum,ne),xm(countnum)
!--out
real(r_size),intent(out)                       :: E(countnum,ne)
!--
integer                                        :: i
!---
do i=1,countnum
    E(:,i)=(x(:,i)-xm(:))
end do
return
end subroutine get_ensemble_diff
!************************************************************************************
! subroutine get_Hxfm()
! end subroutine get_Hxfm
!************************************************************************************
end module varxf