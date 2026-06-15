      ! f2py -c -m read_patchMS read_patchMS.F90 --fcompiler=gnu95
      ! python -m numpy.f2py -c -m read_patchMS read_patchMS.F90
!***************************************************************************
    subroutine read_patchms(ix,iy,threshold,nextX,nextY,nextdst,uparea,&
        & rivseq,nx,ny,xlist,ylist) !,k) outdir,
    !------------------------------------------------
    implicit none
    integer,intent(IN)                    :: ix,iy
    integer,intent(IN)                    :: nx, ny
    real,intent(IN)                       :: threshold
    integer,dimension(nx,ny),intent(IN)   :: nextX,nextY,rivseq
    real,dimension(nx,ny),intent(IN)      :: nextdst,uparea
    !character*128,intent(IN)              :: outdir
    integer,dimension(100),intent(OUT)    :: xlist, ylist
    !integer,intent(OUT)                   :: k

    character(len=8)                      :: llon,llat
    character(len=128)                    :: fname
    real,dimension(nx,ny)                 :: weightage
    integer                               :: info, iix, iiy, fn, k, ios 
    character*128                         :: outdir
    !===========
    ! open emperical weightage
    write(llon,'(i4.4)') ix
    write(llat,'(i4.4)') iy
    outdir='/cluster/data6/menaka/Empirical_LocalPatch'
    fname=trim(adjustl(outdir))//"/weightage/"//trim(llon)//trim(llat)//".bin"
    fn = 34
    open(fn,file=fname,form="unformatted",access="direct",recl=4*ny*nx,status="old",iostat=ios)
    if(ios==0)then
        read(fn,rec=1) weightage
    else
        write(*,*) "no weightage", fname
    end if
    close(fn)
    ! get the mainstream pixels
    ! most upstream pixels and number of downstream pixels
    iix=0
    iiy=0
    k=0
    xlist=-9999
    ylist=-9999
    call patch_pixels(ix,iy,nx,ny,threshold,weightage,nextX,nextY,nextdst,uparea,rivseq,iix,iiy,k)
    !write(*,*)"===",ix,iy,"===",iix,iiy,k
    !allocate(xlist(k),ylist(k))
    call downstream_pixels(iix,iiy,k,nx,ny,nextX,nextY,xlist,ylist,info)
    !print* ,xlist(1:5)
    !print* ,ylist(1:5)
    if (info /= 0 .and. k==0) then
        write(*,*) "info:",info
        xlist(1)=ix
        ylist(1)=iy
        !deallocate(xlist,ylist)
        !cycle
    end if
    return
    end subroutine read_patchMS
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
    !*****************************************************************
    subroutine wgt_consistancy(i,j,x,y,nx,ny,nextX,nextY,weightage,thersold,conflag)
    implicit none 
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: weightage
    !--
    real                                :: conflag,thersold
    integer                             :: ix,iy,iix,iiy,tx,ty,ud
    real                                :: length!,rl
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
      conflag=1.0
      do while (ix/=tx .or. iy/=ty) 
        iix=ix
        iiy=iy 
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        if (ix==-9 .or. iy==-9) then
          ud=+1
          exit
        end if
        if (ix==-10 .or. iy==-10) then
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=+1
          exit
        end if
        !--compare the weightage and thersold
        if (weightage(ix,iy) < thersold) then
          conflag=0.0
        end if 
      end do
    end if
    !--
    if (ud==+1) then
      tx=i
      ty=j
      ix=x
      iy=y
      length=0.0
      conflag=1.0
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
        if (ix==-10 .or. iy==-10) then
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=-9999
          exit
        end if
        !--compare the weightage and thersold
        if (weightage(ix,iy) < thersold) then
          conflag=0.0
        end if 
      end do
    end if
    !-- 
    if (ud==-9999) then
      conflag=0.0
    elseif (ud==0) then
      conflag=1.0
    else
      conflag=conflag
    end if
    !---
    return
    !---
    end subroutine wgt_consistancy
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
    subroutine upstream(i,j,nx,ny,nextX,nextY,uparea,x,y)
    implicit none 
    ! find the upstream pixel with closest upstream area to the i,j
    !--
    integer,intent(IN)                        :: i,j,nx,ny
    integer,dimension(nx,ny),intent(IN)       :: nextX,nextY !,rivseq
    real,dimension(nx,ny),intent(IN)          :: uparea
    integer,intent(OUT)                       ::x,y
    !--
    real                                      :: dA ! area differnce nextdst,
    integer                                   :: ix,iy,iix,iiy,tx,ty,ud,d
    real                                      :: length,rl
    !--
    x=-9999
    y=-9999
    d=100 ! look at 100*100 box
    dA=1.0e20 ! area differnce
    !--
    !write(*,*)i,j
    do tx=i-d,i+d
      do ty=j-d,j+d
        !write(*,*)tx,ty
        call ixy2iixy(tx,ty, nx, ny, ix, iy)
        !write(*,*)nextX(ix,iy),nextY(ix,iy),ix,iy,uparea(ix,iy),rivseq(ix,iy)
        if (nextX(ix,iy) == i .and. nextY(ix,iy) == j) then
            !write(*,*)ix,iy
            if (abs(uparea(i,j)-uparea(ix,iy)) < dA) then
                dA=abs(uparea(i,j)-uparea(ix,iy))
                x=ix
                y=iy
                !write(*,*)x,y
            end if
        end if
      end do
    end do
    return
    end subroutine upstream
    !******************************************************************
    subroutine patch_pixels(i,j,nx,ny,threshold,weight,nextX,nextY,nextdst,uparea,rivseq,x,y,k)
    implicit none
    ! find the upstream pixel with closest upstream area to the i,j
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY,rivseq
    real,dimension(nx,ny)               :: weight,nextdst,uparea
    !--
    real                                :: threshold,dA ! area differnce
    integer                             :: ix,iy,iix,iiy,tx,ty,ud,d,k
    real                                :: length,rl
    !integer,dimension(1000)             :: lx,ly
    !--
    k=0
    !write(*,*)i,j
    x=i
    y=j
    if (rivseq(i,j)>1) then ! no upstream
        ix=i
        iy=j
        do while (weight(ix,iy)>=threshold)
            call upstream(ix,iy,nx,ny,nextX,nextY,uparea,iix,iiy) !rivseq,
            !write(*,*)iix,iiy,rivseq(ix,iy)
            if (weight(iix,iiy)<threshold) exit
            !write(*,*)iix,iiy,weight(iix,iiy),weight(iix,iiy)>=threshold
            if (rivseq(iix,iiy)==1) then
                x=iix
                y=iiy
                exit
            end if
            ix=iix
            iy=iiy
            k=k+1
        end do
        x=ix
        y=iy
    end if
    !write(*,*)"upstream:",k!,weight(x,y)
    ud=k
    k=0
    call downstream(x,y,nx,ny,threshold,weight,nextX,nextY,k)
    !write(*,*)"downstream:",k-ud !,max(k-ud,0)
    return
    end subroutine patch_pixels
    !*****************************************************************
    subroutine downstream(i,j,nx,ny,threshold,weight,nextX,nextY,k)
    implicit none
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: weight !nextdst
    !--
    real                                :: threshold!,lag_dist
    integer                             :: ix,iy,iix,iiy,k!,tx,ty,ud,k
    !real                                :: length,rl
    !--
    k=1
    !--
    ix=i
    iy=j
    do while (weight(ix,iy)>=threshold)
      iix=ix
      iiy=iy 
      ix=nextX(iix,iiy)
      iy=nextY(iix,iiy)
      !write(*,*)ix,iy
      if (ix==-9 .or. iy==-9) then
        k=k+1
        exit
      end if
      if (ix==-10 .or. iy==-10) then
        k=k+1
        exit
      end if
      if (ix==-9999 .or. iy==-9999) then
        k=k+1
        exit
      end if
      k=k+1
    end do
    !---
    k=k-1
    if (k==0) then
        if (weight(i,j)>=threshold) then
            k=1
        end if
    end if
    !write(*,*)k,weight(i,j)
    !---
    return
    !---
    end subroutine downstream
    !*****************************************************************
    subroutine downstream_pixels(i,j,k,nx,ny,nextX,nextY,lx,ly,info)
    implicit none
    !--
    integer                             :: i,j,k,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    !real,dimension(nx,ny)               :: weight !nextdst
    !--
    !real                                :: threshold,lag_dist
    integer                             :: num,iix,iiy,ix,iy,info
    !eal                                :: length,rl
    integer,dimension(k)                :: lx,ly
    !--
    num=1
    ix=i
    iy=j
    lx(num)=ix
    ly(num)=iy
    num=num+1
    !--
    do while (num<=k)
      iix=ix
      iiy=iy
      ix=nextX(iix,iiy)
      iy=nextY(iix,iiy)
      !write(*,*)ix,iy
      if (ix==-9 .or. iy==-9) then
        num=num+1
        exit
      end if
      if (ix==-10 .or. iy==-10) then
        num=num+1
        exit
      end if
      if (ix==-9999 .or. iy==-9999) then
        num=num+1
        exit
      end if
      lx(num)=ix
      ly(num)=iy
      num=num+1
    end do
    !---
    num=num-1
    info=num-k
    !---
    return
    !---
    end subroutine downstream_pixels
    !*****************************************************************