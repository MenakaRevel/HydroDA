program data_assim
!$ use omp_lib
!**************************
! Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019)]
! created by Ikeshima & Menaka
! Menaka@IIS 2020
!**************************
implicit none
character*128                   :: fname,buf,camadir,expdir,DAdir,patchdir
character*8                     :: yyyymmdd,nxtyyyymmdd,befyyyymmdd
real                            :: assimN,assimS,assimW,assimE,lat,lon
character*2                     :: swot_day
character*4                     :: patchid
real,allocatable                :: swot_obs(:,:),global_xa(:,:,:)
integer*4                       :: lon_cent,lat_cent,patch_size,patch_side,i,j,k,countnum,patch_nums,countR
!integer*4                       :: S_lon_cent,S_lat_cent
integer,allocatable             :: local_obs(:),iwork(:),ifail(:),H(:,:)!,localx(:,:,:)
real,allocatable                :: xf_m(:),xf(:,:),globalx(:,:,:),xa(:,:)!,H(:,:)!xa_m(:),,localx_line(:)
real,allocatable                :: meanglobalx(:,:,:),meanglobaltrue(:,:)
integer                         :: ens_num,num,ios,ovs,info,info2,errflg,m
character*3                     :: numch
real,allocatable                :: globaltrue(:,:),xt(:),R(:,:),Rdiag(:)!
real,allocatable                :: W(:,:),Pa(:,:),Pasqr(:,:),UNI(:,:),EfW(:,:),HETRHE(:,:)!,HPH(:,:)
real                            :: VDVTmax!,xf_m55,xa_m55,xt_m55,xf_err,xa_err,delta,traceR,traceHPH
real,allocatable                :: work(:),la(:),U(:,:),Dinv(:,:),VDVT(:,:),Dsqr(:,:),yo(:),Ef(:,:),la_p(:),U_p(:,:)
integer,allocatable             :: isuppz(:)
integer                         :: day!,E_retry
real,allocatable                :: randlist(:)
real                            :: errexp
real,allocatable                :: K_(:,:)

integer,allocatable             :: obs_mask(:,:) ! NEW v.1.1.0
real,allocatable                :: ens_xa(:,:,:)

integer,parameter               :: latpx=720,lonpx=1440
real,dimension(lonpx,latpx)     :: rivwth,rivlen,nextdst,weightage,storage,parm_infl
integer,dimension(lonpx,latpx)  :: nextX,nextY,ocean,countp,targetp

integer,allocatable             :: usedwhat(:)
real                            :: errrand,errfix
!-for inflation
real                            :: rho,rho_fixed ! covariance inflation parameter, 1.01
real                            :: thresold ! weightage thresold , 0.2
real,allocatable                :: dep(:),HEf(:,:),HEfR(:,:),HPH(:,:)
real,dimension(4)               :: parm
real                            :: sigma_o,gain,sigma_b ! background variance
real,parameter                  :: rho_min=1.0d0

!real,allocatable                :: storage(:,:)
real,allocatable                :: global_null(:,:)!,globals_count(:,:)
real,allocatable                :: Wvec(:),lag(:),local_lag(:)!global_sum_xam(:,:),
real,allocatable                :: wgt(:),local_wgt(:)
real                            :: wt,lag_dist!lag_distance,Gauss_wt,
! local variables
real,allocatable                :: local_swot(:),local_err(:)!,local_swot_line(:)
integer,allocatable             :: local_ocean(:)
real,allocatable                :: local_river(:)!,localRW_line(:)
integer*4                       :: i_m,j_m
integer,allocatable             :: xlist(:),ylist(:)
integer*4                       :: target_pixel,fn
character*8                     :: llon,llat
real,dimension(lonpx,latpx)     :: obs_err, obserrrand

write(*,*) "data_assim"
call getarg(1,buf)
read(buf,*) assimN

call getarg(2,buf)
read(buf,*) assimS

call getarg(3,buf)
read(buf,*) assimW

call getarg(4,buf)
read(buf,*) assimE

call getarg(5,buf)
read(buf,*) yyyymmdd

call getarg(6,buf)
read(buf,*) swot_day

call getarg(7,buf)
read(buf,*) patch_size ! radius

call getarg(8,buf)
read(buf,*) ens_num ! number of ensemble

call getarg(9,buf)
read(buf,*) day ! number of date; start from 0

call getarg(10,buf)
read(buf,*) nxtyyyymmdd

call getarg(11,buf)
read(buf,*) befyyyymmdd

call getarg(12,buf)
read(buf,*) errexp

call getarg(13,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(14,buf)
read(buf,*) errrand
write(*,*) errrand

call getarg(15,buf)
read(buf,*) errfix

call getarg(16,buf)
read(buf,*) thresold

call getarg(17,buf)
read(buf,"(A)") expdir

call getarg(18,buf)
read(buf,"(A)") DAdir

call getarg(19,buf)
read(buf,"(A)") patchdir

call getarg(20,buf)
read(buf,*) rho_fixed

call getarg(21,buf)
read(buf,*) sigma_b


allocate(usedwhat(3))
usedwhat=0

write(*,*) errfix

patch_side=patch_size*2+1
patch_nums=patch_side**2

! rho covariance inflation parameter
rho=1.0d0
! weightgae thersold
!thresold=0.8d0
!---
fname=trim(adjustl(expdir))//"/logout/errrand_"//yyyymmdd//".log"
open(36,file=fname,status='replace')
write(36,*) errrand
close(36)

fname=trim(adjustl(expdir))//"/logout/assimLog_"//yyyymmdd//".log"
open(78,file=fname,status='replace')

fname=trim(adjustl(expdir))//"/logout/KLog_"//yyyymmdd//".log"
open(84,file=fname,status='replace')

20 format(i4.4,2x,i4.4,2x,f8.4,2x,f8.4,2x,f8.4)
21 format(i4.4,2x,i4.4,2x,f12.7,2x,f12.7,2x,f12.7)
22 format(a4,2x,a4,2x,a8,2x,a8,2x,a8)
23 format(i4.4,2x,i4.4,2x,f10.7)

fname=trim(adjustl(expdir))//"/logout/testLog"//yyyymmdd//".log"
open(72,file=fname,status='replace')
write(72,22)"lon","lat","true","forcast","assim"

fname=trim(adjustl(expdir))//"/logout/pixelLog_"//yyyymmdd//".log"
open(79,file=fname,status='replace')
write(79,*) "lat","lon","valid pixels in emperical patch"
write(*,*) "lat","lon","valid pixels in emperical patch"

!$ write(*,*)"omp threads",omp_get_num_threads()
fname=trim(adjustl(expdir))//"/logout/inflation_"//yyyymmdd//".log"
open(73,file=fname,status='replace')

fname=trim(adjustl(expdir))//"/logout/ensembles_"//yyyymmdd//".log"
open(74,file=fname,status='replace')

fname=trim(adjustl(expdir))//"/logout/error_"//yyyymmdd//".log"
open(82,file=fname,status='replace')

! read storage (for making ocean mask)
!allocate(ocean(lonpx,latpx))
fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//"T000/storge"//yyyymmdd(1:4)//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) storage
else
    write(*,*) "no file storage"
end if
close(34)

write(*,*) "read storage"

! make ocean mask from storage data (1 is ocean; 0 is not ocean)
!ocean = (storage>1e18) * (-1)

! read river width
fname=trim(adjustl(camadir))//"map/glb_15min/rivwth_gwdlr.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
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
    write(*,*) "no file nextXY at:",fname
end if
close(34)

! make ocean mask from nextx (1 is ocean; 0 is not ocean)
ocean = (nextX==-9999) * (-1)

! read river length
fname=trim(adjustl(camadir))//"map/glb_15min/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)

! read distance to next grid
fname=trim(adjustl(camadir))//"map/glb_15min/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
else
    write(*,*) "no file nextdst",fname
end if
close(34)
!--
! read SWOT observation distance data
allocate(swot_obs(1440,640))
swot_obs=0
fname=trim(adjustl(DAdir))//"/sat/mesh_day"//swot_day//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*640,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) swot_obs
else
    write(*,*) "no swot"
end if
close(34)

! obs_mask is a mask for considering if there is observation in local patch at that date
allocate(obs_mask(lonpx,latpx))
!obs_mask=0 ! 0 means no observations in local patch
!fname=trim(adjustl(DAdir))//"/ava_obs/obs"//swot_day//".bin"
!open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
!if(ios==0)then
!    read(34,rec=1) obs_mask
!else
!    write(*,*) "no obs_mask"
!end if
!close(34)

! inflation parameter
fname=trim(adjustl(expdir))//"/inflation/parm_infl"//yyyymmdd//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) parm_infl
else
    write(*,*) "no parm_infl"
end if
close(34)

! read mean sfelv forcast
allocate(meanglobalx(lonpx,latpx,ens_num))
meanglobalx=0
do num=1,ens_num
    write(numch,'(i3.3)') num
    !fname=trim(adjustl(expdir))//"/assim_out/mean_sfcelv/meansfcelvC"//numch//".bin"
    fname=trim(adjustl(expdir))//"/assim_out/ens_xa/"//befyyyymmdd//"_"//numch//"_xa.bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) meanglobalx(:,:,num)
    else
        write(*,*) "no mean x"
    end if
    close(34)
end do

! read mean WSE true
allocate(meanglobaltrue(lonpx,latpx))
meanglobaltrue=0
!fname=trim(adjustl(expdir))//"/assim_out/mean_sfcelv/meansfcelvT000.bin"
fname=trim(adjustl(expdir))//"/assim_out/mean_sfcelv/meansfcelvT000.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) meanglobaltrue
else
    write(*,*) "no true"
end if
close(34)

! read WSE from all model
allocate(globalx(lonpx,latpx,ens_num))
globalx=0
do num=1,ens_num
    write(numch,'(i3.3)') num
    fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//"A"//numch//"/sfcelv"//yyyymmdd(1:4)//".bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) globalx(:,:,num)
    else
        write(*,*) "no x"
    end if
    close(34)
end do

! update globalx
globalx=globalx-meanglobalx

! mean of ensembles
!do num=1,ens_num
!    globalx(:,:,num)=globalx(:,:,num)-(sum(meanglobalx(:,:,:),dim=3)/real(ens_num))
!end do

! read true WSE
allocate(globaltrue(lonpx,latpx))
globaltrue=0
fname=trim(adjustl(expdir))//"/CaMa_out/"//yyyymmdd//"T000/sfcelv"//yyyymmdd(1:4)//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) globaltrue
else
    write(*,*) "no true"
end if
close(34)
! update globalture
globaltrue=globaltrue-meanglobaltrue!+(sum(meanglobalx(:,:,:),dim=3)/real(ens_num))

! read countnum
!fname=trim(adjustl(expdir))//"/local_patch/countnum.bin"
!fname="../covariance/local_patch_0.90/countnum.bin"
fname=trim(adjustl(patchdir))//"/countnum.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) countp
    read(34,rec=2) targetp
else
    write(*,*) "countp , tsrget"
end if
close(34)

! make global model average
!allocate(global_ave(lonpx,latpx))
!global_ave=0
!do i=1,1440
!    do j=1,720
!        global_ave(i,j)=sum(globalx(i,j,:))/(1e-20+real(ens_num))
!    end do
!end do

! make randomlist
allocate(randlist(ens_num*366))
randlist=0
!fname=trim(adjustl(expdir))//"/CaMa_out/randlist.bin"
!open(34,file=fname,form="unformatted",access="direct",recl=4*366*ens_num,status="old",iostat=ios)
!if(ios==0)then
!    read(34,rec=1) randlist
!else
!    write(*,*) "no randlist"
!end if
!close(34)

! observation error 0.1*(1/L)*(1/W)
fname=trim(adjustl(DAdir))//"/sat/obs_err.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) obs_err
else
    write(*,*) "no obs error", fname
end if
close(34)

! observation error random value 
fname=trim(adjustl(expdir))//"/CaMa_out/errrand.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) obserrrand
else
    write(*,*) "no obs error rand", fname
end if
close(34)

allocate(global_xa(lonpx,latpx,ens_num),global_null(lonpx,latpx))!,global_count(lonpx,latpx)
global_xa = 0
!global_count = 0
global_null = 0.0

! add ERROR to Observation
globaltrue=globaltrue+errrand
!globaltrue=globaltrue+obserrrand

! obs_mask is a mask for considering if there is observation at that date
!allocate(obs_mask(lonpx,latpx))
!obs_mask=0 ! 0 means no observatioin

!write(72,*) "L211"


! parallel calculation

!$omp parallel default(none)&
!$omp& shared(lon_cent,assimW,assimE,assimN,assimS,ocean,rivwth,obs_mask,patch_size,ens_num,&
!$omp& patch_nums,nextX,nextY,nextdst,errfix, &
!$omp& swot_obs,globalx,globaltrue,global_xa,global_null) &
!$omp& private(lat_cent,lat,lon,llon,llat,fname,fn,ios,weightage, &
!$omp& lag,xlist,ylist,wgt,countnum,j,i,i_m,j_m,lag_dist,target_pixel,xt, &
!$omp& local_ocean,local_river,local_lag,local_wgt,local_swot,local_obs, &
!$omp& xf,errflg,ovs,H,Ef,xf_m,R,Rdiag,countR,wt,W,VDVT,Pa,Pasqr, &
!$omp& UNI,la_p,U_p,HETRHE,VDVTmax,work,iwork,ifail,m,info,U,la,Dinv,Dsqr,info2,yo, &
!$omp& Wvec,xa,EfW,K_,num,rho)
!$omp do
do lon_cent = int((assimW+180)*4+1),int((assimE+180)*4+1),1
  do lat_cent = int((90-assimN)*4+1),int((90-assimS)*4+1),1
        lat = 90.0-(lat_cent-1.0)/4.0
        lon = (lon_cent-1.0)/4.0-180.0
        !write(*,*)"================",lat,lon!,omp_get_num_threads(),"================="
        ! not connected to longtitude direction; no calculation available near lon=-180,180 or lat=-80,80
        !remove ocean
        if (ocean(lon_cent,lat_cent)==1) then
            cycle
            !continue
        end if
        ! remove rivwth <= 0m
        if (rivwth(lon_cent,lat_cent) <=0.0) then
            cycle
        end if
        ! observation available in local patch
!        if (obs_mask(lon_cent,lat_cent)==0) then
!            cycle
!        end if
        !write(*,*)"===",lat,lon,"==="!$,omp_get_num_threads(),omp_get_thread_num(),"==="
        !write(78,*)lat,lon
        ! make xt =======================================
        !write(*,*) "allocate lag"
        !allocate(lag(patch_nums))!localRW_line(patch_nums))
        !write(*,*) "allocate xlist ylist"
        ! countnum and target pixel
        countnum=countp(lon_cent,lat_cent)
        target_pixel=targetp(lon_cent,lat_cent)

!25 format(countnum*(i4.4,2x,i4.4,2x,f10.7))

        !write(*,*)countnum, target_pixel
        ! open emperical local patch
        allocate(lag(patch_nums),xlist(countnum),ylist(countnum),wgt(countnum))
        write(llon,'(i4.4)') lon_cent
        write(llat,'(i4.4)') lat_cent
        !fname="./local_patch/patch"//trim(llon)//trim(llat)//".txt"
        !fname="../covariance/local_patch_0.90/patch"//trim(llon)//trim(llat)//".txt"
        fname=trim(adjustl(patchdir))//"/patch"//trim(llon)//trim(llat)//".txt"
        !write(*,*) fname
        open(34,file=fname,status='old',access='sequential',form='formatted',action='read')!
        do i=1, countnum
          read(34,*) xlist(i),ylist(i),wgt(i)
          !write(*,*) xlist(i),ylist(i),wgt(i)
        end do
        !write(*,23) xlist,ylist,wgt
        close(34)
        !--
        allocate(xt(countnum),local_ocean(countnum),local_river(countnum),local_lag(countnum),local_wgt(countnum))
        !--
        !local_lag=lag(1:countnum)
        local_wgt=wgt(1:countnum)
        !--
        !write(*,*) xlist !"make xt"
        ! creating local ocean and river pixels
        xt=0
        local_ocean=0
        local_river=0.0
        do i=1,countnum
            i_m=xlist(i)
            j_m=ylist(i)
            xt(i)=globaltrue(i_m,j_m)
            local_ocean(i)=ocean(i_m,j_m)
            local_river(i)=rivwth(i_m,j_m)
        end do
        !--
        ! read local SWOT distances
        allocate(local_swot(countnum),local_err(countnum))
        local_swot=0.0
        local_err=0.0
        !--
        do i=1,countnum
            i_m=xlist(i)
            j_m=ylist(i)
            local_err(i)=obs_err(i_m,j_m)
            if (j_m<=40) then
                local_swot(i)=-9999.
            elseif (j_m>680) then
                local_swot(i)=-9999.
            else
                local_swot(i)=swot_obs(i_m,j_m-40)
            end if
        end do

        !local_swot_line=local_swot!reshape(local_swot,(/patch_nums/))
        !deallocate(local_swot)

        allocate(local_obs(countnum))
        local_obs=0

        ! SWOT
        local_obs=(local_swot<=60. .and. local_swot>=10.)*(-1) ! .true.=1 of .false.=0
        !write(72,*) lat,lon,local_obs

        ! ocean is excluded
        local_obs=local_obs*(local_ocean==0)*(-1)

        ! river width<=50 is excluded
        local_obs=local_obs*(local_river>50)*(-1)


        ! make xf =====================================
        !write(*,*) "make xf"
        allocate(xf(countnum,ens_num))!localx(countnum,countnum,ens_num),
        !allocate(local_river(patch_side,patch_side),localRW_line(patch_nums))
        !localx=0
        !localx_line=0
        xf=0
        !write(*,*) "make xf"
        !do k=1,ens_num ! k: ensemble number
        do i=1,countnum
            i_m=xlist(i)
            j_m=ylist(i)
            xf(i,:)=globalx(i_m,j_m,:)
        end do
            !xf(:,k)=localx_line(:)
        !end do
        ! make xf_m ====================================
        !write(*,*) "make xf_m"

        ! deallocate variables of making SWOT observation and dimension related
        ! variables
        deallocate(lag,local_ocean,local_river,local_swot,xlist,ylist,wgt)

        errflg=0
        ! calculate the number of observations
        if(sum(local_obs)==0)then
            !xa=xf
            errflg=1
            !write(*,*) "error",errflg
            write(82,*) lat,lon,"error",errflg
            goto 9999
        end if
        !------------------------------------------------------------
        ! ========= reach here only when there is observation =======
        !------------------------------------------------------------
        write(78,*) "================================================"
        !write(78,*) "******************",lat,lon,"*******************"
        write(78,*) "******************",lon_cent,lat_cent,"*******************"
        write(*,*) "******************",lon_cent,lat_cent,"*******************"
        write(78,*) "size",countnum
        write(78,*) "local obs",local_obs

        ! inflation parameter
        if (rho_fixed==-1.0) then
            rho=parm_infl(lon_cent,lat_cent)
        else
            rho=rho_fixed
        endif
        if (rho<rho_min) rho=rho_min
        write(*,*)rho

        ! number observations
        ovs=sum(local_obs)
        if(ovs>0)then
            ! observation available
            allocate(H(ovs,countnum))
            H=0
            j=1 ! row number
            do i=1,countnum ! col number
                if(local_obs(i)==1)then
                    H(j,i)=1
                    j=j+1
                end if
            end do
        end if

        !write(78,*) "obs available: ",lat,lon

        ! count how many pixels were valid in the local patch
        !write(79,*) lat,lon,sum(local_obs)
        !write(78,*) lat,lon,sum(local_obs)

        ! make obs_mask 1
        !obs_mask(lon_cent,lat_cent)=1

        ! make Ef ======================================
        allocate(Ef(countnum,ens_num),xf_m(countnum))
        !allocate(xf_m(countnum))
        Ef=0
        xf_m=0

        do i=1,countnum
            xf_m(i)=sum(xf(i,:))/(1e-20+real(ens_num))
        end do

        do k=1,ens_num
            Ef(:,k)=(xf(:,k)-xf_m(:))
        end do

        write(78,*) "shape H:",shape(H)
        !write(78,*) "H:",H
        !write(78,*) "Ef:",Ef
        !write(78,*) "xf:",xf
        !write(78,*) "xf_m:",xf_m
        !write(78,*) "xt:",xt
        !write(78,*) "wt:", matmul(H,local_wgt)
        !write(78,*) "lag_dist",local_lag*local_obs

        write(79,*) "================================================"
        write(79,*) lon_cent,lat_cent !lat,lon
        write(79,*) "wt:", matmul(H,local_wgt)

        ! make R (NEW: add weigtage) ============================
        ! added observation localization
        allocate(R(ovs,ovs),Rdiag(ovs))
        countR=1
        !errfix=obs_err(lon_cent,lat_cent)
        !if (errfix > 0.25) then ! assume largest uncertinity is 0.25 m
        !    errfix=0.25
        !end if
        do i=1, countnum
          if (local_obs(i)==1) then
              !wt=Gauss_wt(local_lag(i))
              !errfix=local_err(i)
              wt=local_wgt(i)
              !wt=max(1.0,wt*1.5)
              Rdiag(countR)=(errfix**2.0)/wt
              countR=countR+1
          end if
        end do

        !write(78,*) "obs:",ovs
        !write(78,*) "Rdiag:",Rdiag

        R=0.
        do i=1,ovs
          R(i,i)=(Rdiag(i))**(-1.)
        end do

        !if(ovs>0) then
        !    write(78,*) "R:",(1/(R+1e-20))
        !    write(78,*) "Ef:",(sum(abs(Ef))/(float(ens_num)+1e-20))**2.
        !end if

        ! make W ====================================
        allocate(W(ens_num,ens_num),VDVT(ens_num,ens_num))
        allocate(Pa(ens_num,ens_num),Pasqr(ens_num,ens_num),UNI(ens_num,ens_num))
        allocate(la_p(ens_num),U_p(ens_num,ens_num),HETRHE(ens_num,ens_num))
        W=0
        VDVT=0
        Pa=0
        Pasqr=0
        UNI=0
        la_p=0
        U_p=0

        UNI=RESHAPE([(1,(0,i=1,ens_num),j=1,ens_num-1),1],[ens_num,ens_num])
        HETRHE=matmul(matmul(TRANSPOSE(matmul(H,Ef)),R),matmul(H,Ef))
        VDVTmax=maxval(abs(HETRHE))

        !write(78,*) "HETRHE:",HETRHE
        !write(78,*) "ET*E:",matmul(TRANSPOSE(Ef),Ef)

        allocate(work(1000),iwork(1000),ifail(1000),isuppz(1000))
        work=0
        iwork=0
        ifail=0
        info=0
        isuppz=0

        ! calculate VDVT
        !if (rho_fixed/=-1.0) then
        !    VDVT=real(ens_num-1.)*UNI/rho_fixed+HETRHE
        !else
        !    VDVT=real(ens_num-1.)*UNI/rho+HETRHE
        !end if
        VDVT=real(ens_num-1.)*UNI/rho+HETRHE
        !   ! modify VDVT
        !   do i=1,ens_num
        !       do j=1,ens_num
        !           if(i<j) VDVT(i,j)=0.
        !       end do
        !   end do

        ! Eigon value decomposition
        ! Diagonalizate VDVT to calculate inverse matrix
        call ssyevx("V","A","U",ens_num,VDVT,ens_num,1e-5,1e5,1,2,1.2e-38*2.,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)

        !call ssyevr("V","A","U",ens_num,VDVT,ens_num,1e-5,1e5,1,2,1.2e-38*2.,m,la_p,U_p,ens_num,isuppz,work,1000,iwork,1000,info)
        !write(78,*) "m",m
        !write(78,*) "ovs:",ovs
        !write(78,*) "la_p",la_p
        !write(78,*) "U_p",U_p
        !write(78,*) "info",info


        if (m<ens_num)then
            errflg=2
            !write(78,*) "~~~ m<ens_num ~~~ m=",m
            write(*,*) "~~~ m<ens_num ~~~ m=",m,rho
            !xa=xf
            write(36,*) lat,lon,"error",errflg,"info",info
            goto 9999 
        end if
        allocate(U(m,m),la(m))
        U=0
        la=0
        do i=1,m
            do j=1,m
                U(i,j)=U_p(i,j)
            end do
            la(i)=la_p(i)
        end do
        allocate(Dinv(m,m),Dsqr(m,m))
        Dinv=0
        Dsqr=0

        ! calc Dinv,Dsqr
        if(info==0)then
            do i=1,m
                Dinv(i,i)=(la(i)+1e-20)**(-1)
            end do
            !=write(78,*) "Dinv",Dinv
            Dsqr=Dinv
            info2=0
            call spotrf("U",m,Dsqr,m,info2)
            !if(info2/=0) write(78,*) "====== ERROR cannot unpack Dsqr======"
            if(info2/=0) then
                errflg=3
                !xa=xf
                write(*,*) "====== ERROR cannot unpack Dsqr======",errflg
                write(82,*) lat,lon,"error",errflg
                goto 9999
            end if
            !write(78,*) "Dsqr",Dsqr
            !write(78,*) "Dinv",Dinv
            Pa   =matmul(matmul(U,Dinv),transpose(U))
            Pasqr=matmul(matmul(U,Dsqr),transpose(U))

            allocate(yo(ovs),Wvec(ens_num))
            yo=matmul(H,xt)

            Wvec = matmul(matmul(matmul(Pa,TRANSPOSE(matmul(H,Ef))),R),yo-matmul(H,xf_m))
            !write(78,*) "yo-Hxt",yo-matmul(H,xf_m)

            do i = 1,ens_num
                W(:,i) = Wvec + sqrt(ens_num-1.)*Pasqr(:,i)
            end do

            deallocate(Wvec)
        else
            !write(78,*) "NG INFO"
            write(*,*) "NG INFO"
            W=0
            errflg=4
            write(82,*) lat,lon,"error",errflg,"info:",info
            goto 9999
        end if


        ! make xa ====================================
        allocate(xa(countnum,ens_num))
        xa=0
        allocate(EfW(countnum,ens_num))
        EfW=0

        EfW=matmul(Ef,W)
        !=write(78,*) "EfW",EfW
        do i=1,ens_num
            xa(:,i)=EfW(:,i)+xf_m(:)
        end do

        !write(79,*) yo-matmul(H,xf_m)
        !write(78,*) "xa:",xa
        !write(*,*) "xa:",xa
        ! check center pixel ====================================
        write(*,*) "errfix:", errfix, obserrrand(lon_cent,lat_cent)
        write(*,*) "true   :",xt(target_pixel)
        write(*,*) "forcast:",sum(xf(target_pixel,:))/(ens_num+1e-20)
        write(*,*) "assimil:",sum(xa(target_pixel,:))/(ens_num+1e-20)


        write(78,*) "true   :",xt(target_pixel)
        write(78,*) "forcast:",sum(xf(target_pixel,:))/(ens_num+1e-20)
        write(78,*) "assimil:",sum(xa(target_pixel,:))/(ens_num+1e-20)
        ! check K_ value (should be between 0-1) =======================================
        allocate(K_(countnum,ovs))
        K_ = matmul(Ef,matmul(matmul(Pa,TRANSPOSE(matmul(H,Ef))),R))
        !write(78,*) "K:",K_
        write(78,*) "K:",sum(K_)
        !write(84,*) "K:",K_
        !write(72,21) lon_cent,lat_cent,xt(target_pixel), sum(xf(target_pixel,:))/(ens_num+1e-20),sum(xa(target_pixel,:))/(ens_num+1e-20)
        write(72,*) lon_cent,lat_cent,xt(target_pixel), sum(xf(target_pixel,:))/(ens_num+1e-20),sum(xa(target_pixel,:))/(ens_num+1e-20)
        write(74,*) "+++++++++++++++++++++++++++++++++++++"
        write(74,*) lon_cent,lat_cent
        write(74,*) "true   :", xt(target_pixel)
        do i=1, ens_num
             write(74,*) "ensemble",i,"forcast:",xf(target_pixel,i),"assimil:",xa(target_pixel,i)
        end do
        !=====================================
        !-- cal inflation parameter estimation
        parm=0.0d0
        allocate(dep(ovs),HEf(ovs,ens_num),HEfR(ens_num,ovs),HPH(ovs,ovs))
        !write(*,*) yo-matmul(H,xf_m)
        dep=yo-matmul(H,xf_m)
        HEf=matmul(H,Ef)
        HEfR=matmul(TRANSPOSE(HEf),R)
        HPH=matmul(HEf,TRANSPOSE(HEf))
        do i=1,ovs
          parm(1)=parm(1)+dep(i)*dep(i)*Rdiag(i)**(-1.)
        end do
        !do num=1,ens_num
        do i=1,ovs
          parm(2)=parm(2)+HPH(i,i)*Rdiag(i)**(-1.)
        end do
        !end do
        parm(2)=parm(2)/real(ens_num-1)
        parm(3)=sum(local_wgt)
        parm(4)=(parm(1)-parm(3))/(parm(2)+1.0e-20) - rho
        sigma_o=2.0d0/(parm(3)+1.0e-20) * ((rho*parm(2) + parm(3))/(parm(2)+1.0e-20))**2
        gain=sigma_b**2 / (sigma_o + sigma_b**2)
        write(73,*) "+++++++++++++++++++++++++++++++++++"
        write(73,*) gain, parm(4),parm(1),parm(2),parm(3)
        rho=rho+ gain* parm(4)
        write(*,*)"rho",rho,rho_min
        if (rho<rho_min) rho=rho_min
        write(73,*) lon_cent,lat_cent,"rho",rho,"rho_min",rho_min
        !---
        parm_infl(lon_cent,lat_cent)=rho
        !write(*,*) target_pixel,shape(xa), xa(target_pixel,num)
        do num=1,ens_num
            global_xa(lon_cent,lat_cent,num) = xa(target_pixel,num)
        end do
        global_null(lon_cent,lat_cent) = 1.0
!        if (sum(K_) > real(ens_num)) then 
!            global_null(lon_cent,lat_cent) = 0.0
!        end if
!
!        if (sum(K_) < 0.0) then 
!            global_null(lon_cent,lat_cent) = 0.0
!        end if
        9999 continue


        ! clean memory ====================================
        ! deallocateとかやる
        if(errflg==0)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,EfW,H,iwork,ifail,U_p,la_p,HETRHE,xf_m,xa,yo,dep,HEf,HEfR,HPH,isuppz)
            !deallocate(local_roff,local_flow,RR_rate)
            deallocate(K_)
        elseif(errflg==2)then
            deallocate(Ef,R,Rdiag,W,VDVT,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m,isuppz)
        elseif(errflg==3)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m,isuppz)
        elseif(errflg==4)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m,isuppz)
        end if
        deallocate(local_obs,xf,xt,local_lag,local_wgt,local_err)
    end do
end do
!$omp end do
!$omp end parallel

fname=trim(adjustl(expdir))//"/logout/usedwhat_"//yyyymmdd//".log"
!open(77,file=fname,status='replace')
!write(77,*) "usedwhat ensemble mean/observation/analysed",usedwhat
!close(77)

! for locations out of assimilation range, make global_null = 1
  ! also, make global_count 1
!   global_null = 0
!   global_null = (global_count<0.5)*(-1)

! calculate global_xam from global_sum_xam and global_count
!   global_xam = global_sum_xam / (global_count+1e-20)

! if global_null == 1
!   global_xam = global_xam*(global_null==0)*(-1) + global_ave*(global_null==1)*(-1)


! NEW v.1.1.0
! make ensemble output (WSE)
allocate(ens_xa(lonpx,latpx,ens_num))
do num=1,ens_num
    ens_xa(:,:,num) = global_xa(:,:,num)*(global_null) + globalx(:,:,num)*(1-global_null) + meanglobalx(:,:,num) !+ (sum(meanglobalx(:,:,:),dim=3)/real(ens_num)) !
end do

!fname="./logout/OutSfcLog_"//yyyymmdd//".log"
!open(85,file=fname,status='replace')
!! output result log
!do i = 1,1440
!    do j = 1,720
!        if(global_null(i,j)==1) write(85,*) &
!        globaltrue(i,j), ens_xa(i,j,1), globalx(i,j,1)
!    end do
!end do
!close(85)

! output ens_xa (assimilated + each ensemble result) NEW v.1.1.0
do num=1,ens_num
    write(numch,'(i3.3)') num
    fname=trim(adjustl(expdir))//"/assim_out/ens_xa/assim/"//yyyymmdd//"_"//numch//"_xa.bin"
    write(*,*)fname
    open(35,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
    if(ios==0)then
        write(35,rec=1) ens_xa(:,:,num)
    else
        write(*,*) "not created", fname
    end if
    close(35)
end do

! save inflation parameter 
fname=trim(adjustl(expdir))//"/inflation/parm_infl"//nxtyyyymmdd//".bin"
open(35,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(35,rec=1) parm_infl
else
    write(*,*) "no parm_infl"
end if
close(35)


!write(*,*) "L.531"

close(78)
close(79)
close(84)
close(72)
close(73)
close(74)
close(82)
deallocate(global_xa,swot_obs,globalx,globaltrue,ens_xa,global_null,obs_mask)
deallocate(meanglobalx,meanglobaltrue)
end program data_assim
!*****************************************************************
subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
real                                :: lag_dist
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length,rl
!--
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
!--
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  lag_dist=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!--
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    !---
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!-- 
if (ud==-9999) then
  lag_dist=-9999
elseif (ud==0) then
  lag_dist=0.0
else
  lag_dist=length
end if
!---
return
!---
end subroutine lag_distance
!**************************************************
subroutine read_wgt(fname,nx,ny,weightage)
!$ use omp_lib    
implicit none
character*128                      :: fname 
integer                            :: fn,nx,ny,ios
real,dimension(nx,ny)              :: weightage
fn=34
!$ fn= fn + omp_get_thread_num()
!!$ write(*,*) fn
open(fn,file=fname,form="unformatted",access="direct",recl=4*ny*nx,status="old",iostat=ios)
if(ios==0)then
    read(fn,rec=1) weightage
else
    write(*,*) "no weightage", fname
end if
close(fn)
!--
return
!---
end subroutine read_wgt
!**************************************************
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
!***************************************************   
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     ix, nx
!-- for output ----------
integer                     roundx
!------------------------
if (ix .ge. 1) then
  roundx = ix - int((ix -1)/nx)*nx
else
  roundx = nx - abs(mod(ix,nx))
end if 
return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
implicit none
!- for input -----------------
integer                   ix, iy, nx, ny
!- for output ----------------
integer                   iix, iiy,roundx
!-----------------------------
if (iy .lt. 1) then
  iiy = 2 - iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else if (iy .gt. ny) then
  iiy = 2*ny -iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else
  iiy = iy
  iix = roundx(ix, nx)
end if
return
end subroutine ixy2iixy
!*****************************************************************        
