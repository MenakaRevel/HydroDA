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
real,intent(in)                                :: globalx(nx,ny,ne)!,meanglobalx(nx,ny,ne),stdglobalx(nx,ny,ne)
!--out
!integer,intent(out)                            :: local_sat(countnum)
real,intent(out)                               :: xf(countnum,ne) !,local_err(countnum)
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
subroutine get_ensemble_mean(x,countnum,ne,Xm)
!=======================================================================
! get ensemble mean of x
!=======================================================================
!--in
integer,intent(in)                             :: ne,countnum
real,intent(in)                                :: X(countnum,ne)
!--out
real,intent(out)                               :: Xm(countnum)
!--
integer                                        :: i
!---
do i=1,countnum
    Xm(i)=sum(X(i,:))/(1e-20+real(ne))
end do
return
end subroutine get_ensemble_mean
!************************************************************************************
subroutine get_ensemble_diff(X,Xm,countnum,ne,E)
!=======================================================================
! get ensemble mean differnce of x
!=======================================================================
!--in
integer,intent(in)                             :: ne,countnum
real,intent(in)                                :: X(countnum,ne),Xm(countnum)
!--out
real,intent(out)                               :: E(countnum,ne)
!--
integer                                        :: i
!---
do i=1,countnum
    E(:,i)=(x(:,i)-xm(:))
end do
return
end subroutine get_ensemble_diff
!************************************************************************************
subroutine get_HX(globalvar,local_obs,xlist,ylist,Hobs,nx,ny,nvar,nobs,ne,countnum,HX)
!=======================================================================
! get HX - simulations in obervational space with ensembles
! input  
!   globalvar - global simulated variable: globalvar[nx,ny,numvar,ne]
!   local_obs - array of available observation: local_obs[nvar*countnum]; 1=observation available
!   xlist     - list of x corrdinates
!   ylist     - list of y corrdinates
!   Hobs      - convert observation to only observation space: Hobs[nobs,countnum]
!   nx        - x dimension of global map
!   ny        - y dimension of global map
!   nvar      - number of variables
!   nobs      - total numer of observations: sum(wse,dis,wsa)
!   ne        - number of ensembles
!   countnum  - number pixel in the local patch
! output
!   HX        - simulations in obervational space: HX[nobs,ne]
!=======================================================================
!--in
integer,intent(in)                             :: nx,ny,numvar,nobs,ne
real,intent(in)                                :: xlist(countnum),ylist(countnum),local_obs(nvar*countnum)
real,intent(in)                                :: globalvar(nx,ny,nvar,ne),Hobs(nobs,countnum)
!--out
real,intent(out)                               :: HX(nobs,ne)
!--
integer                                        :: i,i_m,j_m,var
! real                                           :: Xt
HX=0
j=1
do var=1, numvar
    do i=1, countnum
        i_m=xlist(i)
        j_m=ylist(i)
        if (local_obs(j) == 1.0) then
            HX(j,:)=globalvar(i_m,j_m,var,:)
        end if
    end do
end do
end subroutine get_HX
!************************************************************************************
end module varxf