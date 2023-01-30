module obser
!====================================================================================
! purpose: Procedures realted to observations: reading observations, error, etc
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
subroutine read_observation(yyyymmdd,nx,ny,obs,obs_err,mean_obs,std_obs)
implicit none
!-in
character(len=8),intent(in)         :: yyyymmdd
integer,intent(in)                  :: nx,ny
!--out
real,dimension(nx,ny),intent(out)   :: obs,obs_err,mean_obs,std_obs
!--
integer                             :: ix,iy,ios
character(len=128)                  :: fname,sat
real                                :: wse,mean,std,obs_error
!--
obs=-9999.0
obs_err=-9999.0
mean_obs=-9999.0
std_obs=-9999.0
    fname="./assim_out/obs/"//trim(yyyymmdd)//".txt"
    print*, fname
    open(11, file=fname, form='formatted',iostat=ios)
    if (ios /= 0) then 
        print*, "no observations: ", fname, ios
        goto 1090
    end if
1000 continue
    read(11,*,end=1090) ix, iy, wse, mean, std, sat
    ! print*, yyyymmdd, ix, iy, wse, trim(sat)
    obs(ix,iy)=wse
    obs_err(ix,iy)=observation_error(sat)
    mean_obs(ix,iy)=mean
    std_obs(ix,iy)=std
    goto 1000
1090 continue
return
end subroutine read_observation
!************************************************************************************
function str2int(str) result(int)
implicit none
! for input
character(len=*),intent(in)            :: str
! for output
integer                                :: int
!-
integer                                :: stat
!---
read(str,*,iostat=stat)  int
if (stat/=0) int=-9999
end function str2int
!************************************************************************************
function observation_error(sat) result(obs_error)
implicit none
!---in
character(len=128),intent(in)           :: sat
!---out
real                                    :: obs_error
!==================================================
! observation errors are from Breada et al,. (2019)
! Brêda, J. P. L. F., Paiva, R. C. D., Bravo, J. M., Passaia, O. A., & Moreira, D. M. (2019). 
! Assimilation of Satellite Altimetry Data for Effective River Bathymetry. Water Resources Research, 
! 55(9), 7441–7463. doi:10.1029/2018WR024010
!--------------------------------------------------
! Satellite  |  Error (m)
! -----------------------
!ENVISAT     |   0.30
!JASON2      |   0.28
!JASON3      |   0.28
!SENTINEL3A  |   0.30*
!SWOT        |   0.10+
!------------------------
! *1.Watson, C.; Legresy, B.; King, M.; Deane, A. Absolute altimeter bias results from Bass Strait, 
! Australia.In Proceedings of the Ocean Surface Topography Science Team Meeting 2018, Ponta Delgada, 
! Portugal,24–29 September 2018.
! *2.Bonnefond, P.; Exertier, P.; Laurain, O.; Guinle, T.; Féménias, P. Corsica:  
! A 20-Yr multi-mission absolute altimeter calibration site.  In Proceedings of the Ocean Surface 
! Topography Science Team Meeting 2018,Ponta Delgada, Portugal, 24–29 September 2018.
! *3.Garcia-Mondejar, A.; Zhao, Z.; Rhines, P. Sentinel-3 Range and Datation Calibration with Crete transponder.
! In Proceedings of the 25 Years of Progress in Radar Altimetry, Ponta Delgada, Portugal, 24–29 September 2018.
! *4.Mertikas, S.; Donlon, C.; Féménias, P.; Mavrocordatos, C.; Galanakis, D.; Tripolitsiotis, A.; Frantzis, X.; 
! Kokolakis, C.; Tziavos, I.N.; Vergos, G.; Guinle, T. Absolute Calibration of the European Sentinel-3A Surface 
! Topography Mission over the Permanent Facility for Altimetry Calibration in west Crete, Greece. Remote Sens. 2018, 10, 1808.
! =================================================
obs_error=0.30
if (trim(sat) == "ENVISAT") obs_error=0.30
if (trim(sat) == "JASON2") obs_error=0.28
if (trim(sat) == "JASON3") obs_error=0.28
if (trim(sat) == "SENTINEL3A") obs_error=0.30
if (trim(sat) == "SENTINEL3B") obs_error=0.30
if (trim(sat) == "SWOT") obs_error=0.10
end function observation_error
!************************************************************************************
subroutine read_local_obs(xlist,ylist,conflag,obs,obs_err,mean_obs,std_obs,countnum,patch_start,patch_end,nx,ny,local_sat,xt,local_err)
!=======================================================================
! get local observations 
! convert observations to differnt space
!!! observation converstions 
!  1 - Directly values 
!  2 - Anomalies
!  3 - Normalized values
!  4 - Log converted values
!=======================================================================
implicit none
!--in
integer,intent(in)                             :: countnum,patch_start,patch_end,nx,ny,conflag
integer,intent(in)                             :: xlist(countnum),ylist(countnum)
real,intent(in)                                :: obs(nx,ny),obs_err(nx,ny),mean_obs(nx,ny),std_obs(nx,ny)
!--out
integer,intent(out)                            :: local_sat(countnum)
real,intent(out)                               :: xt(countnum),local_err(countnum)
!--
integer                                        :: i,j,i_m,j_m
!============
! initialize
local_sat=-9999.0
xt=-9999.0
local_err=-9999.0
j=1
do i=patch_start,patch_end
    i_m=xlist(i)
    j_m=ylist(i)
    if (obs(i_m,j_m)/=-9999.0) then
        local_sat(j)=1.0
        ! convert observations
        if (conflag == 1) then
            xt(j)=obs(i_m,j_m)
            local_err(j)=obs_err(i_m,j_m)
        else if (conflag == 2) then
            xt(j)=obs(i_m,j_m)-mean_obs(i_m,j_m)
            local_err(j)=obs_err(i_m,j_m)
        else if (conflag == 3) then
            xt(j)=(obs(i_m,j_m)-mean_obs(i_m,j_m))/(std_obs(i_m,j_m)+1.0e-20)
            local_err(j)=obs_err(i_m,j_m)/(std_obs(i_m,j_m)+1.0e-20)
        else if (conflag == 4) then
            xt(j)=log10(obs(i_m,j_m))
            local_err(j)=sqrt(log10(obs_err(i_m,j_m)**2+1))
        end if
    end if
    j=j+1
end do
return
end subroutine read_local_obs
!************************************************************************************
end module obser