program calc_stat
implicit none
character*128                   :: fname,buf,camadir,outdir,assim_out
character*8                     :: yyyymmdd
character*3                     :: numch
integer*4                       :: lon_cent,lat_cent
integer,parameter               :: latpx=720,lonpx=1440
integer,dimension(lonpx,latpx)  :: nextX,nextY,ocean
real,dimension(lonpx,latpx)     :: AI,RMSE,NRMSEasm,NRMSEopn,rRMSE,VEasm,VEopn,NSA,NSC,NSE,PDRI,PTRI,ENSPR,PBIASasm,PBIASopn,KGEasm,KGEopn ! ensemble spread
real,allocatable                :: opn(:,:,:,:),asm(:,:,:,:),org(:,:,:)
real,allocatable                :: opn_mean(:),asm_mean(:),org_mean(:),err_mask(:)
integer                         :: ios,N,m,i,day
integer,allocatable             :: days(:)
real                            :: dis_mean
real                            :: Qtp, Qap, Qcp
integer                         :: Ttp,Tap,Tcp,Tp0,Tp1
real,allocatable                :: opn_max(:),opn_min(:),asm_max(:),asm_min(:)
real                            :: KGE

call getarg(1,buf)
read(buf,*) N ! number of days in year 366/365
write(*,*)N

call getarg(2,buf)
read(buf,*) m  !number of ensembles

call getarg(3,buf)
read(buf,"(A)") outdir
write(*,*) outdir

call getarg(4,buf)
read(buf,"(A)") camadir
write(*,*) camadir

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/glb_15min/nextxy.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextX
    read(34,rec=2) nextY
else
    write(*,*) "no file nextXY at:",fname
end if
close(34)
!-----
write(*,*)"ocean"
ocean=(nextx==-9999)*(-1)
!-----
! read names of days in a particular year
allocate(days(N))
fname="year_day.txt"
open(34,file=fname,form="formatted",iostat=ios)
if(ios==0)then
    read(34,*) days
else
    write(*,*) "no days",fname
end if
close(34)
!--
allocate(org(N,lonpx,latpx),opn(N,m,lonpx,latpx),asm(N,m,lonpx,latpx))
org=0.0
opn=0.0
asm=0.0
write(*,*)"read rivout"
do day=1,N
    write(yyyymmdd,'(i8.0)') days(day)
    write(*,*) yyyymmdd !day,, days(day)
    ! true river dischrge
    fname=trim(adjustl(outdir))//"/assim_out/rivout/true/rivout"//yyyymmdd//".bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) org(day,:,:)
    else
        write(*,*) "no true discharge",fname
    end if
    close(34)
    do i=1,m
        write(numch,'(i3.3)') i
        ! corrpted river discharge
        fname=trim(adjustl(outdir))//"/assim_out/rivout/open/rivout"//yyyymmdd//"_"//numch//".bin"
        open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
        if(ios==0)then
            read(34,rec=1) opn(day,i,:,:)
        else
            write(*,*) "no corrupted discharge",fname
        end if
        close(34)
        ! assimilated river discharge
        fname=trim(adjustl(outdir))//"/assim_out/rivout/assim/rivout"//yyyymmdd//"_"//numch//".bin"
        open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
        if(ios==0)then
            read(34,rec=1) asm(day,i,:,:)
        else
            write(*,*) "no assimilated discharge",fname
        end if
        close(34)
    end do
end do
!--
allocate(org_mean(N),opn_mean(N),asm_mean(N),err_mask(N))
allocate(opn_max(N),opn_min(N),asm_max(N),asm_min(N))
do lon_cent = 1,lonpx
    do lat_cent = 1, latpx
        !remove ocean
        if (ocean(lon_cent,lat_cent)==1) then
            cycle
            !continue
        end if
        ! AI calculation
        org_mean=org(:,lon_cent,lat_cent)
        opn_mean=sum(opn(:,:,lon_cent,lat_cent),dim=2)/real(m)
        asm_mean=sum(asm(:,:,lon_cent,lat_cent),dim=2)/real(m)
        !write(*,*) shape(asm_mean) , N
        err_mask=((abs(org_mean-asm_mean)/(org_mean+1.0e-20))>0.1)*(-1)
        !write(*,*)sum(err_mask)
        !write(*,*)sum(org_mean)/real(N)
        AI(lon_cent,lat_cent)= sum((1.0 -(abs(org_mean-asm_mean)/abs(org_mean-opn_mean+1.0e-20))),mask=err_mask==0.0)/(sum(err_mask)+1.0e-20)
        ! RMSE
        RMSE(lon_cent,lat_cent)=sqrt((1/real(N))*sum((asm_mean-org_mean)**2))
        ! rRMSE
        rRMSE(lon_cent,lat_cent)=sqrt((1/real(N))*sum(((asm_mean-org_mean)/(org_mean+1.0e-20))**2))
        ! NRMSE
        dis_mean=sum(org_mean)/(real(N))
        NRMSEasm(lon_cent,lat_cent)=sqrt((1/real(N))*sum((asm_mean-org_mean)**2))/(dis_mean+1.0e-20)
        NRMSEopn(lon_cent,lat_cent)=sqrt((1/real(N))*sum((opn_mean-org_mean)**2))/(dis_mean+1.0e-20)
        ! VEasm
        VEasm(lon_cent,lat_cent)=1.0 - (sum(asm_mean-org_mean)/(sum(org_mean)+1.0e-20))
        ! VEopn
        VEopn(lon_cent,lat_cent)=1.0 - (sum(opn_mean-org_mean)/(sum(org_mean)+1.0e-20))
        ! NSA
        NSA(lon_cent,lat_cent)=1.0 - (sum((asm_mean-org_mean)**2)/(sum((dis_mean-org_mean)**2)+1.0e-20))
        ! NSC
        NSC(lon_cent,lat_cent)=1.0 - (sum((opn_mean-org_mean)**2)/(sum((dis_mean-org_mean)**2)+1.0e-20))
        ! NSE
        NSE(lon_cent,lat_cent)=((NSA(lon_cent,lat_cent)-NSC(lon_cent,lat_cent))/(1-NSC(lon_cent,lat_cent)))
        if (NSE(lon_cent,lat_cent) == 0.0) then
            if ( NSA(lon_cent,lat_cent)/= NSC(lon_cent,lat_cent) ) then
                write(*,*)NSA(lon_cent,lat_cent),NSC(lon_cent,lat_cent),(NSA(lon_cent,lat_cent)-NSC(lon_cent,lat_cent)),1-NSC(lon_cent,lat_cent)
            end if
        end if
            ! PDRI & PTRI
        ! peak discharge and peak timing
        Qtp=maxval(org_mean)
        Ttp=maxloc(org_mean,dim=1)
        ! peak discharge and timing of assimilated and corrupted
        Tp0=max(Ttp-15,0)
        Tp1=min(Ttp+15,N)
        Qap=maxval(asm_mean(Tp0:Tp1))
        Tap=maxloc(asm_mean(Tp0:Tp1),dim=1)
        Qcp=maxval(opn_mean(Tp0:Tp1))
        Tcp=maxloc(opn_mean(Tp0:Tp1),dim=1)
        ! PDRI
        PDRI(lon_cent,lat_cent)=1.0 - (abs(Qtp-Qap)/(abs(Qtp-Qcp)+1.0e-20))
        if (Qtp==Qcp) then
            PDRI(lon_cent,lat_cent)=0.0
        end if
        ! PTRI
        PTRI(lon_cent,lat_cent)=1.0 - (real(abs(Ttp-Tap))/(real(abs(Ttp-Tcp))+1.0e-20))
        if (Ttp==Tcp) then
            PTRI(lon_cent,lat_cent)=0.0
        end if
        ! Ensemble Spread
        opn_max=maxval(opn(:,:,lon_cent,lat_cent),dim=2)
        opn_min=minval(opn(:,:,lon_cent,lat_cent),dim=2)
        asm_max=maxval(asm(:,:,lon_cent,lat_cent),dim=2)
        asm_min=minval(asm(:,:,lon_cent,lat_cent),dim=2)
        !
        ENSPR(lon_cent,lat_cent)=sum(1.0 -((asm_max-asm_min)/(opn_max-opn_min+1.0e-20)))/(real(N))
        !
        PBIASasm(lon_cent,lat_cent)=100.0*(sum(asm_mean-org_mean)/(sum(org_mean)+1.0e-20))
        PBIASopn(lon_cent,lat_cent)=100.0*(sum(opn_mean-org_mean)/(sum(org_mean)+1.0e-20))
        ! KGE
        KGEasm(lon_cent,lat_cent)=KGE(asm_mean,org_mean,N)
        KGEopn(lon_cent,lat_cent)=KGE(opn_mean,org_mean,N)
        !write(*,*)lon_cent,lat_cent,AI(lon_cent,lat_cent),NSA(lon_cent,lat_cent),NSC(lon_cent,lat_cent)!((NSA(lon_cent,lat_cent)-NSC(lon_cent,lat_cent))/(1-NSC(lon_cent,lat_cent))),KGEasm(lon_cent,lat_cent)
       !,rRMSE(lon_cent,lat_cent),NRMSEasm(lon_cent,lat_cent),VE(lon_cent,lat_cent),NSE(lon_cent,lat_cent),,PBIASasm(lon_cent,lat_cent)
    end do
end do
deallocate(org,opn,asm)
deallocate(org_mean,opn_mean,asm_mean,err_mask)
deallocate(days)
deallocate(opn_max,opn_min,asm_max,asm_min)

! remove AI<0.0
!AI=AI*((AI>0.0)*(1.0))

! NSE
!NSE= ((NSA-NSC)/(1.0-NSC))
! remove NSE<0.0
!NSE=NSE*((NSE>0.0)*(1.0))


! Assimilation Index
fname=trim(adjustl(outdir))//"/assim_out/stat/annualmeanAI.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) AI
else
    write(*,*) "no AI",fname
end if
close(34)

! Ensemble Spread
fname=trim(adjustl(outdir))//"/assim_out/stat/annualmeanENSPR.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1)ENSPR
else
    write(*,*) "no Ensemble Spread",fname
end if
close(34)

! RMSE
fname=trim(adjustl(outdir))//"/assim_out/stat/RMSEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) RMSE
else
    write(*,*) "no RMSE",fname
end if
close(34)

! rRMSE
fname=trim(adjustl(outdir))//"/assim_out/stat/rRMSEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) rRMSE
else
    write(*,*) "no rRMSE",fname
end if
close(34)

! NRMSE
fname=trim(adjustl(outdir))//"/assim_out/stat/NRMSEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) NRMSEasm
else
    write(*,*) "no NRMSEasm",fname
end if
close(34)

fname=trim(adjustl(outdir))//"/assim_out/stat/NRMSEopn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) NRMSEopn
else
    write(*,*) "no NRMSEopn",fname
end if
close(34)

! VEasm
fname=trim(adjustl(outdir))//"/assim_out/stat/VEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) VEasm
else
    write(*,*) "no VE assimilated",fname
end if
close(34)

! VE
fname=trim(adjustl(outdir))//"/assim_out/stat/VEopn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) VEopn
else
    write(*,*) "no VE corrupted",fname
end if
close(34)


! NSEasm
fname=trim(adjustl(outdir))//"/assim_out/stat/NSEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) NSA
else
    write(*,*) "no NSA",fname
end if
close(34)

! NSEopn
fname=trim(adjustl(outdir))//"/assim_out/stat/NSEopn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) NSC
else
    write(*,*) "no NSC",fname
end if
close(34)

! NSEAI
fname=trim(adjustl(outdir))//"/assim_out/stat/NSEAI.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) NSE
else
    write(*,*) "no NSE",fname
end if
close(34)

! PDRI
fname=trim(adjustl(outdir))//"/assim_out/stat/annualPDRI.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) PDRI
else
    write(*,*) "no PDRI",fname
end if
close(34)

! PTRI
fname=trim(adjustl(outdir))//"/assim_out/stat/annualPTRI.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) PTRI
else
    write(*,*) "no PDRI",fname
end if
close(34)

! pBias assimilated
fname=trim(adjustl(outdir))//"/assim_out/stat/pBIASasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) PBIASasm
else
    write(*,*) "no pBias assimilated",fname
end if
close(34)

! pBias corrupted
fname=trim(adjustl(outdir))//"/assim_out/stat/pBIASopn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) PBIASopn
else
    write(*,*) "no pBias corrupted",fname
end if
close(34)

! KGE corrupted
fname=trim(adjustl(outdir))//"/assim_out/stat/KGEopn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) KGEopn
else
    write(*,*) "no KGE corrupted",fname
end if
close(34)

! KGE assimilated
fname=trim(adjustl(outdir))//"/assim_out/stat/KGEasm.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(34,rec=1) KGEasm
else
    write(*,*) "no KGE simulated",fname
end if
close(34)

end program calc_stat
!**********************
function std(list,N)
implicit none
!--
integer               :: N
real,dimension(N)     :: list
!--
integer               :: i
real                  :: var,std,mean
!--
var=0.0d0
mean=sum(list)/(real(N)+1.0e-20)
!--
do i=1, N
    var=var+(list(i)-mean)**2
end do
!--
std=sqrt((var)/(real(N)+1.0e-20))
!--
return
end function
!************************
function cov(X,Y,N)
implicit none
!--
integer                :: N
real,dimension(N)      :: X,Y
real                   :: cov
!--
real                   :: X_mean,Y_mean
real                   :: C
integer                :: i
!--
X_mean=sum(X)/(real(N)+1.0e-20)
Y_mean=sum(Y)/(real(N)+1.0e-20)
!--
C=0.0d0
!--
do i=1,N
    C=C+(X(i)-X_mean)*(Y(i)-Y_mean)
end do
cov=C/((N-1)+1.0e-20)
!--
return
!--
end function cov
!************************
function KGE(sim,obs,N)
implicit none
!--
integer                :: N
real,dimension(N)      :: sim,obs
real                   :: KGE
!--
real                   :: sim_mean,obs_mean,sim_std,obs_std
real                   :: std,cov
real                   :: CC,BR,RV
!--
sim_mean=sum(sim)/(real(N)+1.0e-20)
obs_mean=sum(obs)/(real(N)+1.0e-20)
!--
sim_std=std(sim,N)
obs_std=std(obs,N)
!--
CC=cov(sim,obs,N)/((sim_std*obs_std)+1.0e-20)
BR=sim_mean/(obs_mean+1.0e-20)
RV=((sim_std/sim_mean)/((obs_std/obs_mean)+1.0e-20))
!--
KGE=1-sqrt((CC-1)**2+(BR-1)**2+(RV-1)**2)
return
!--
end function KGE
!*************************
