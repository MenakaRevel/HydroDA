program data_assim
!$ use omp_lib
! ====================================================================================
! ADD LICENSE
! ====================================================================================
!*************************************************************************************
! Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019,2021)]
! ====================================================================================
! Reference:
! 1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating 
! global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 
! 1–34. https://doi.org/10.1029/2020wr027876
! 2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical 
! Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: 
! A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
! ====================================================================================
! version 1 was created by Ikeshima & Menaka
! Menaka@IIS 2021
!
! revised and modified by Menaka
! Menaka@MSU 2026
!*************************************************************************************
    use letkf
    use obser
    use patch
    use varxf
    use common
    use prase
    !*************************************************************************************
    implicit none

    ! --- File I/O and String Variables ---
    character(len=128)              :: fname, buf, camadir, expdir, DAdir, patchdir, hydrowebdir, mapname, patchname, cal
    character(len=8)                :: yyyymmdd, nxtyyyymmdd
    character(len=3)                :: numch
    character(len=256)              :: cvarsin, varsda, string_buffer
    character(len=256)              :: varnames(100), obstypes(100)
    character(len=64)               :: cvnames_in(100), cvnames_da(100)
    character(len=128)              :: station            ! VS name
    character(len=4)                :: llon, llat
    integer                         :: file_unit, io_status, ios, ovs, info, info2, errflg, m
    
    ! --- Domain and Map Variables ---
    real(r_size)                    :: assimN, assimS, assimW, assimE, lat, lon
    real(r_size)                    :: gsize, west, north, east, south ! map boundaries
    integer                         :: ny, nx, nflp              ! pixel size, calculated
    real(r_size), allocatable       :: rivwth(:,:), rivhgt(:,:), rivlen(:,:), nextdst(:,:), lons(:,:), lats(:,:)
    real(r_size), allocatable       :: elevtn(:,:), fldhgt(:,:,:)
    integer, allocatable            :: nextX(:,:), nextY(:,:), ocean(:,:), countp(:,:), targetp(:,:)
    
    ! --- Local Patch Variables ---
    integer(kind=4)                 :: lon_cent, lat_cent, patch_size, patch_side, patch_nums, countR
    integer(kind=4)                 :: patch_start, patch_end, countnumber, targetpixel, target_pixel, countnum
    integer(kind=4)                 :: i, j, k, i_m, j_m
    real(r_size), allocatable       :: Wvec(:), lag(:), local_lag(:)
    real(r_size), allocatable       :: wgt(:), local_wgt(:)
    real(r_size)                    :: wt, lag_dist
    integer, allocatable            :: xlist(:), ylist(:)
    
    ! --- Observation Variables ---
    integer                         :: nobstypes, nobsvars, nda, nvar, nvars, flag
    real(r_size)                    :: wse, std            ! observed wse and std
    real(r_size), allocatable       :: obs(:,:), obs_err(:,:), altitude(:,:), mean_obs(:,:), std_obs(:,:)
    real(r_size), allocatable       :: local_sat(:), local_err(:), swot_obs(:)
    integer, allocatable            :: local_obs(:), vobs(:)
    integer                         :: conflag, nobs, obs_mask
    integer, allocatable            :: conflags(:)
    real(r_size)                    :: pslamch
    
    ! --- Assimilation & State Variables ---
    integer                         :: ens_num, num, day, ne
    integer                         :: lwork, liwork ! added for large ensemble size
    integer                         :: local_ocean, local_river
    real(r_size)                    :: sigma_o, gain, sigma_b ! background variance
    real(r_size)                    :: infl, min_infl, infl_flg
    real(r_size), parameter         :: rho_min = 1.0d0
    real(r_size)                    :: rho, rho_fixed ! covariance inflation parameter
    real(r_size)                    :: thresold       ! weightage thresold
    real(r_size)                    :: errrand, errfix, errexp, VDVTmax
    real(r_size), allocatable       :: global_xa(:,:,:), global_null(:,:) 
    real(r_size), allocatable       :: globalx(:,:,:), globalhxb(:,:,:,:), ens_xa(:,:,:)
    real(r_size), allocatable       :: meanglobalx(:,:,:), stdglobalx(:,:,:), meanglobaltrue(:,:), stdglobaltrue(:,:)
    real(r_size), allocatable       :: xf_m(:), xf(:,:), xa(:,:), globaltrue(:,:), xt(:), R(:,:)
    real(r_size), allocatable       :: Rdiag(:), Rwgt(:), T(:,:)
    real(r_size), allocatable       :: W(:,:), Pa(:,:), Pasqr(:,:), UNI(:,:), EfW(:,:), HETRHE(:,:)
    real(r_size), allocatable       :: work(:), la(:), U(:,:), Dinv(:,:), VDVT(:,:), Dsqr(:,:), Ef(:,:), la_p(:), U_p(:,:)
    real(r_size), allocatable       :: yo(:), dep(:), HEfR(:,:), HPH(:,:), HEf(:,:), K_(:,:)
    real(r_size), allocatable       :: Hobs(:,:), HXb(:,:)
    real(r_size), allocatable       :: weightage(:,:), storage(:,:), parm_infl(:,:)
    real(r_size), dimension(4)      :: parm
    integer, allocatable            :: iwork(:), ifail(:), H(:,:), isuppz(:), usedwhat(:)

    ! -------------------------------------------------------------
    ! Group Variables into the Namelist
    ! -------------------------------------------------------------
    namelist /config_vars/ &
        mapname, patch_size, ens_num, &
        camadir, thresold, expdir, DAdir, patchdir, patchname, &
        hydrowebdir, rho_fixed, sigma_b, conflag, cal, &
        varsda, obstypes

    ! -------------------------------------------------------------
    ! Read the namelist
    ! -------------------------------------------------------------
    open(newunit=file_unit, file='input.nml', status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        write(*,*) "Error: Could not open input.nml file!"
        stop
    end if
    read(file_unit, nml=config_vars, iostat=io_status)
    if (io_status /= 0) then
        write(*,*) "Error: Failed to read namelist config_vars!"
        stop
    end if
    close(file_unit)

    write(*,*) "data_assim"

    call getarg(1, buf)
    read(buf,*) yyyymmdd

    call getarg(2, buf)
    read(buf,*) nxtyyyymmdd

    ! --- Read Map Parameters ---
    fname = trim(camadir) // "/map/" // trim(mapname) // "/params.txt"
    open(11, file=fname, form='formatted')
    read(11,*) nx
    read(11,*) ny
    read(11,*) nflp
    read(11,*) gsize
    read(11,*) west
    read(11,*) east
    read(11,*) south
    read(11,*) north
    close(11)

    ! --- Update the assimilation domain ---
    assimN = min(north,   80.0)
    assimS = max(south,  -60.0)
    assimW = max(west , -180.0)
    assimE = min(east ,  180.0)
    print*, assimN, assimS, assimW, assimE
    
    ! --- Formatting ---
20  format(i4.4,2x,i4.4,2x,f8.4,2x,f8.4,2x,f8.4)
21  format(i4.4,2x,i4.4,2x,f12.7,2x,f12.7,2x,f12.7)
22  format(a4,2x,a4,2x,a8,2x,a8,2x,a8)
23  format(i4.4,2x,i4.4,2x,f10.7)

    allocate(usedwhat(3))
    usedwhat = 0

    print *, "ensemble number", ens_num

    patch_side = patch_size * 2 + 1
    patch_nums = patch_side**2

    rho = 1.0d0 ! rho covariance inflation parameter

    lwork = max(1, 26 * ens_num)
    liwork = max(1, 10 * ens_num)
    
    ! --- Initiate logfiles ---
    fname = trim(adjustl(expdir)) // "/logout/errrand_" // yyyymmdd // ".log"
    open(36, file=fname, status='replace')
    errrand = -1
    write(36,*) errrand
    close(36)
    
    fname = trim(adjustl(expdir)) // "/logout/assimLog_" // yyyymmdd // ".log"
    open(78, file=fname, status='replace')
    write(78,*) "File I/O Errors"

    fname = trim(adjustl(expdir)) // "/logout/KLog_" // yyyymmdd // ".log"
    open(84, file=fname, status='replace')

    fname = trim(adjustl(expdir)) // "/logout/testLog" // yyyymmdd // ".log"
    open(72, file=fname, status='replace')
    write(72,22) "lon", "lat", "true", "forcast", "assim"

    fname = trim(adjustl(expdir)) // "/logout/pixelLog_" // yyyymmdd // ".log"
    open(79, file=fname, status='replace')
    write(79,*) "lat", "lon", "valid pixels in emperical patch"
    write(*,*) "lat ", "lon ", "valid pixels in emperical patch"

    fname = trim(adjustl(expdir)) // "/logout/inflation_" // yyyymmdd // ".log"
    open(73, file=fname, status='replace')

    fname = trim(adjustl(expdir)) // "/logout/ensembles_" // yyyymmdd // ".log"
    open(74, file=fname, status='replace')

    ! I/O related error will be written in /logout/error_{yyyymmdd}.log
    fname = trim(adjustl(expdir)) // "/logout/error_" // yyyymmdd // ".log"
    open(82, file=fname, status='replace')
    write(82,*) "File I/O Errors"

    ! --- Allocate Map Arrays ---
    allocate(rivwth(nx,ny), rivhgt(nx,ny), fldhgt(nx,ny,10), rivlen(nx,ny))
    allocate(nextdst(nx,ny), lons(nx,ny), lats(nx,ny))
    allocate(elevtn(nx,ny), weightage(nx,ny), storage(nx,ny), parm_infl(nx,ny))
    allocate(nextX(nx,ny), nextY(nx,ny), ocean(nx,ny), countp(nx,ny), targetp(nx,ny))

    ! read river width
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/rivwth_gwdlr.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) rivwth
    else
        write(*,*) "no file rivwth"
        write(82,*) "no file rivwth at:", fname
        write(78,*) "no file rivwth at:", fname
    end if
    close(34)

    ! read river channel depth
    if (trim(cal) == "yes") then
        fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/rivhgt_Xudong.bin"
    elseif (trim(cal) == "corrupt") then
        fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/rivhgt_corrupt.bin"
    else
        fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/rivhgt.bin"
    end if
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) rivhgt
    else
        write(*,*) "no file rivhgt"
        write(82,*) "no file rivhgt at:", fname
        write(78,*) "no file rivhgt at:", fname
    end if
    close(34)

    ! read flood plain height
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/fldhgt.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx*10, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) fldhgt
    else
        write(*,*) "no file fldhgt"
        write(82,*) "no file fldhgt at:", fname
        write(78,*) "no file fldhgt at:", fname
    end if
    close(34)

    ! read nextX and nextY
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/nextxy.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) nextX
        read(34, rec=2) nextY
    else
        write(*,*) "no file nextXY at:", fname
        write(82,*) "no file nextXY at:", fname
        write(78,*) "no file nextXY at:", fname
    end if
    close(34)

    ! read lons and lats
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/lonlat.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) lons
        read(34, rec=2) lats
    else
        write(*,*) "no file lonlat at:", fname
        write(82,*) "no file lonlat at:", fname
        write(78,*) "no file lonlat at:", fname
    end if
    close(34)

    ! make ocean mask from nextx (1 is ocean; 0 is not ocean)
    ocean = merge(-1, 0, nextX <= 0)

    ! read river length
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/rivlen.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) rivlen
    else
        write(*,*) "no file rivlen at", fname
        write(82,*) "no file rivlen at:", fname
        write(78,*) "no file rivlen at:", fname
    end if
    close(34)

    ! read distance to next grid
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/nxtdst.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) nextdst
    else
        write(*,*) "no file nextdst at", fname
        write(82,*) "no file nextdst at:", fname
        write(78,*) "no file nextdst at:", fname
    end if
    close(34)

    ! read elevation data
    fname = trim(adjustl(camadir)) // "/map/" // trim(mapname) // "/elevtn.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) elevtn
    else
        write(*,*) "no file elevtn at", fname
        write(82,*) "no file elevtn at:", fname
        write(78,*) "no file elevtn at:", fname
    end if
    close(34)

    ! read observations and observation error variance
    allocate(obs(nx,ny), obs_err(nx,ny), altitude(nx,ny), mean_obs(nx,ny), std_obs(nx,ny))

    ! obstypes
    write(78,*) "========================================================="
    print*, "read observations"
    write(78,*) "read observations"
    call read_observation(yyyymmdd, obstypes, nx, ny, nvars, obs, obs_err, mean_obs, std_obs)

    ! inflation parameter
    fname = trim(adjustl(expdir)) // "/inflation/parm_infl" // yyyymmdd // ".bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*nx*ny, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) parm_infl
    else
        write(*,*) "no parm_infl", fname
        write(82,*) "no file parm_infl:", fname
        write(78,*) "no file parm_infl at:", fname
    end if
    close(34)

    ! read mean sfelv forcast
    allocate(meanglobalx(nx,ny,ens_num), stdglobalx(nx,ny,ens_num))
    meanglobalx = 0
    do num = 1, ens_num
        write(numch, '(i3.3)') num
        fname = trim(adjustl(expdir)) // "/assim_out/mean_sfcelv/meansfcelvC" // numch // ".bin"
        open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
        if(ios == 0) then
            read(34, rec=1) meanglobalx(:,:,num)
        else
            write(*,*) "no mean x", fname
            write(82,*) "no mean x at:", fname
            write(78,*) "no mean x at:", fname
        end if
        close(34)
    end do

    stdglobalx = 0
    do num = 1, ens_num
        write(numch, '(i3.3)') num
        fname = trim(adjustl(expdir)) // "/assim_out/mean_sfcelv/stdsfcelvC" // numch // ".bin"
        open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
        if(ios == 0) then
            read(34, rec=1) stdglobalx(:,:,num)
        else
            write(*,*) "no std x", fname
            write(82,*) "no std x at:", fname
            write(78,*) "no std x at:", fname
        end if
        close(34)
    end do

    ! read mean and std WSE true
    allocate(meanglobaltrue(nx,ny))
    meanglobaltrue = 0.0
    fname = trim(adjustl(expdir)) // "/assim_out/mean_sfcelv/mean_sfcelv.bin"

    allocate(stdglobaltrue(nx,ny))
    stdglobaltrue = 0.0
    fname = trim(adjustl(expdir)) // "/assim_out/mean_sfcelv/std_sfcelv.bin"

    !*************************************************************************************
    ! read water storage - prognostic variable from all model
    !*************************************************************************************
    allocate(globalx(nx,ny,ens_num))
    globalx = 0
    do num = 1, ens_num
        write(numch, '(i3.3)') num
        fname = trim(adjustl(expdir)) // "/CaMa_out/" // yyyymmdd // "A" // numch // "/storge" // yyyymmdd(1:4) // ".bin" 
        open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
        if(ios == 0) then
            read(34, rec=1) globalx(:,:,num)
        else
            write(*,*) "no x :", fname
            write(82,*) "no x at:", fname
            write(78,*) "no x at:", fname
        end if
        close(34)
    end do
    
    !=======================================================================
    ! Conversion flags for variables
    !=======================================================================
    call parse_csv_string(varsda, 250, nvars, varnames)
    allocate(conflags(nvars))
    do nvar = 1, nvars
        if (varnames(nvar) == 'sfcelv') then
            conflags(nvar) = 2 ! anomaly
        elseif (varnames(nvar) == 'outflw') then
            conflags(nvar) = 4 ! lg converted
        elseif (varnames(nvar) == 'fldara') then
            conflags(nvar) = 5 ! precentage innudated
        elseif (varnames(nvar) == 'fldwdh') then 
            conflags(nvar) = 3 ! normalized value
        else
            conflags(nvar) = 1 ! direct value
        end if
    end do

    !=======================================================================
    ! read CMF variables
    !=======================================================================
    allocate(globalhxb(nx,ny,nvars,ens_num)) 
    globalhxb = 0  
    do nvar = 1, nvars
        do num = 1, ens_num
            write(numch, '(i3.3)') num
            fname = trim(adjustl(expdir)) // "/CaMa_out/" // yyyymmdd // "A" // &
                  numch // "/" // trim(varnames(nvar)) // yyyymmdd(1:4) // ".bin"
            if(ios == 0) then
                read(34, rec=1) globalhxb(:,:,nvar,num) 
            else
                write(*,*) "no x :", fname
                write(82,*) "no x at:", fname
                write(78,*) "no x at:", fname
            end if
            close(34)
        end do
    end do

    !=======================================================================
    ! Read the local patch related variables
    !=======================================================================
    fname = trim(adjustl(patchdir)) // "/" // trim(patchname) // "/countnum.bin"
    open(34, file=fname, form="unformatted", access="direct", recl=4*ny*nx, status="old", iostat=ios)
    if(ios == 0) then
        read(34, rec=1) countp
        read(34, rec=2) targetp
    else
        write(*,*) "no file :countp , target", fname
        write(82,*) "no countnum at:", fname
        write(78,*) "no countnum at:", fname
    end if
    close(34)

    !=======================================================================
    ! Allocate global x
    !=======================================================================
    allocate(global_xa(nx,ny,ens_num), global_null(nx,ny))
    global_xa = 0
    global_null = 0.0

    write(82,*) "====================================="
    write(82,*) "Calculation Errors"
    write(78,*) "====================================="
    write(78,*) "Assimilation of each grid"

! parallel calculation
!$omp parallel default(none) &
!$omp shared(assimW, assimE, assimN, assimS, west, east, north, south, gsize, &
!$omp ocean, rivwth, rivhgt, patch_size, ens_num, patch_nums, &
!$omp nextX, nextY, nextdst, globalx, globaltrue, global_xa, global_null, &
!$omp lats, lons, countp, targetp, patchdir, patchname, conflags, &
!$omp obs, obs_err, mean_obs, std_obs, nx, ny, parm_infl) &
!$omp private(lon_cent, lat_cent, lat, lon, llon, llat, fname, ios, weightage, &
!$omp countnumber, targetpixel, lag, xlist, ylist, wgt, &
!$omp patch_start, patch_end, target_pixel, countnum, local_sat, local_err, &
!$omp local_lag, local_wgt, local_obs, xt, nvar, vobs, xf, Hobs, &
!$omp nobs, infl, ne, HEf, Rdiag, Rwgt, yo, HXb, min_infl, infl_flg, T, errflg, &
!$omp EfW, i, xa, xf_m, Ef, j, i_m, j_m, lag_dist, local_ocean, local_river, &
!$omp ovs, H, R, countR, wt, W, VDVT, Pa, Pasqr, UNI, la_p, U_p, HETRHE, &
!$omp VDVTmax, work, iwork, ifail, m, info, U, la, Dinv, Dsqr, info2, Wvec, K_, num, rho)
!$omp do
    do lon_cent = int((assimW-west)*(1.0/gsize)+1), int((assimE-west)*(1.0/gsize)), 1
        do lat_cent = int((north-assimN)*(1.0/gsize)+1), int((north-assimS)*(1.0/gsize)), 1
            lat = lats(lon_cent,lat_cent) 
            lon = lons(lon_cent,lat_cent) 
            
            if (ocean(lon_cent,lat_cent) == 1) cycle
            if (rivwth(lon_cent,lat_cent) <= 0.0) cycle
            
            !=========================
            write(82,*) "+++++++++++++++++"
            write(82,*) lon_cent, lat_cent
            
            countnumber = countp(lon_cent,lat_cent)
            targetpixel = targetp(lon_cent,lat_cent)
            
            if (targetpixel == -9999) cycle
            
            allocate(lag(patch_nums), xlist(countnumber), ylist(countnumber), wgt(countnumber))
            
            write(llon, '(i4.4)') lon_cent
            write(llat, '(i4.4)') lat_cent

            !============================
            ! read emperical local patch 
            !============================
            fname = trim(adjustl(patchdir)) // "/" // trim(patchname) // "/patch" // trim(llon) // trim(llat) // ".txt"
            call read_elp(fname, countnumber, xlist, ylist, wgt)
            
            ! prepare local patch start and end
            call assign_local_patch(countnumber, targetpixel, patch_size, patch_start, patch_end, target_pixel, countnum)
            
            write(79,*) "patch dimesion", patch_start, patch_end, target_pixel, countnum
            
            allocate(local_sat(countnum), local_err(countnum), local_lag(countnum))
            allocate(local_wgt(countnum), local_obs(countnum), xt(countnum))
            
            local_sat = -9999.0
            local_err = -9999.0
            local_wgt = wgt(patch_start:patch_end)
            write(79,*) "local_wgt", local_wgt
            xt = -9999.0

            !====================================================
            ! read local observations
            !====================================================
            call read_local_obs(xlist, ylist, conflags, obs, obs_err, mean_obs, std_obs, countnum, patch_start, patch_end, nx, ny, nvar, local_sat, xt, local_err, vobs)
            
            local_obs = merge(-1, 0, local_sat /= -9999.0)
            
            if (sum(local_obs) == 0) then
                errflg = 1
                goto 9999
            end if

            write(78,*) "=========================================================="
            write(78,*) "******************", lon_cent, lat_cent, " *******************"
            write(78,*) "=========================================================="
            print*, "******************", lon_cent, lat_cent, "*******************"
            write(78,*) "size", countnum
            write(78,*) "Processing..."

            !====================================================
            ! read local prognostic variable
            !====================================================
            allocate(xf(countnum,ens_num))
            xf = 0
            call local_xf(globalx, xlist, ylist, countnum, patch_start, patch_end, nx, ny, ens_num, xf)
            
            call get_Hobs(local_sat, countnum, nvar, vobs, nobs, Hobs)

            !===============================================================
            ! do data assimilation here
            !===============================================================
            infl = parm_infl(lon_cent,lat_cent)
            call letkf_core(ne, nobs, HEf, Rdiag, Rwgt, Yo, HXb, infl, min_infl, infl_flg, T, errflg)

            if (errflg /= 0) then
                write(*,*) "LETKF Core Failed with error code: ", errflg
                xa = xf 
            else
                EfW = matmul(Ef, T)
                do i = 1, ens_num
                    xa(:, i) = xf_m(:) + EfW(:, i)
                end do
                parm_infl(lon_cent,lat_cent) = infl
                write(78,*) "forcast:", sum(xf(target_pixel,:)) / (ens_num + 1e-20)
                write(78,*) "assimil:", sum(xa(target_pixel,:)) / (ens_num + 1e-20)
            end if

9999        continue
            ! Explicit Deallocations to prevent runtime memory errors on next iteration
            if (allocated(lag)) deallocate(lag)
            if (allocated(xlist)) deallocate(xlist)
            if (allocated(ylist)) deallocate(ylist)
            if (allocated(wgt)) deallocate(wgt)
            if (allocated(local_sat)) deallocate(local_sat)
            if (allocated(local_err)) deallocate(local_err)
            if (allocated(local_lag)) deallocate(local_lag)
            if (allocated(local_wgt)) deallocate(local_wgt)
            if (allocated(local_obs)) deallocate(local_obs)
            if (allocated(xt)) deallocate(xt)
            if (allocated(xf)) deallocate(xf)
        end do
    end do
!$omp end do
!$omp end parallel

end program data_assim