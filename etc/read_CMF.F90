! f2py -c -m read_CMF read_CMF.F90 --fcompiler=gnu95
! python -m numpy.f2py -c -m read_CMF read_CMF.F90
! created by Xudong Zhou @IIS
! ************************************************
subroutine read_discharge(ix1, iy1, ix2, iy2, ndays, foutflw,nXX, nYY, dis)
! ===============================================
implicit none
integer,intent(IN)                :: ix1,iy1
integer,intent(IN)                :: ix2,iy2
integer,intent(IN)                :: ndays
integer,intent(IN)                :: nXX, nYY
character*128,intent(IN)              :: foutflw
real,dimension(ndays),intent(OUT) :: dis

real,dimension(nXX,nYY)           :: outflw
integer                           :: irec


open(12, file=foutflw, form='unformatted', access='direct', &
    & recl=4*nXX*nYY)

do irec = 1,ndays
    read(12,rec=irec) outflw
    if (ix2 .GT. 0 .and. iy2 .GT. 0 ) then
    dis(irec)=outflw(ix1,iy1) + outflw(ix2,iy2)
    else
    dis(irec)=outflw(ix1,iy1)
    end if 
    !print *, irec, dis(irec)
end do

close(12)
end subroutine read_discharge
! ************************************************
subroutine read_discharge_multi(sno,ix1, iy1, ix2, iy2, ndays, foutflw,nXX, nYY, dis)
! ===============================================
implicit none
integer,intent(IN)                                :: sno
integer,dimension(sno),intent(IN)                :: ix1,iy1
integer,dimension(sno),intent(IN)                :: ix2,iy2
integer,intent(IN)                :: ndays
integer,intent(IN)                :: nXX, nYY
character*128,intent(IN)              :: foutflw
real,dimension(sno,ndays),intent(OUT) :: dis

real,dimension(nXX,nYY)           :: outflw
integer                           :: irec
integer                           :: isno, ix1_sno,ix2_sno,iy1_sno,iy2_sno


open(12, file=foutflw, form='unformatted', access='direct', &
    & recl=4*nXX*nYY)

do irec = 1,ndays
    read(12,rec=irec) outflw
    do isno = 1, sno
        ix1_sno = ix1(isno)+1
        ix2_sno = ix2(isno)+1
        iy1_sno = iy1(isno)+1
        iy2_sno = iy2(isno)+1
        if (ix2_sno .GT. 0 .and. iy2_sno .GT. 0 ) then
        dis(isno,irec)=outflw(ix1_sno,iy1_sno) + outflw(ix2_sno,iy2_sno)
        else
        dis(isno,irec)=outflw(ix1_sno,iy1_sno)
        end if 
        !print *, irec, dis(irec)
    end do
enddo

close(12)
end subroutine read_discharge_multi
! ************************************************
subroutine read_sfcelv(ix, iy, ndays, fsfcelv,nXX, nYY, wse)
! ===============================================
        implicit none
        integer,intent(IN)                    :: ix,iy
        !integer,intent(IN)                   :: ix2,iy2
        integer,intent(IN)                    :: ndays
        integer,intent(IN)                    :: nXX, nYY
        character*128,intent(IN)              :: fsfcelv
        real,dimension(ndays),intent(OUT)     :: wse

        real,dimension(nXX,nYY)               :: sfcelv
        integer                               :: irec


        open(12, file=fsfcelv, form='unformatted', access='direct', &
            & recl=4*nXX*nYY)
        
        do irec = 1,ndays
            read(12,rec=irec) sfcelv
            if (ix .GT. 0 .and. iy .GT. 0 ) then
            wse(irec)=sfcelv(ix,iy)
            ! else
            !   wse(irec)=outflw(ix1,iy1)
            end if 
            !print *, irec, wse(irec)
        end do

        close(12)
        end subroutine read_sfcelv
!***************************************************************************
        subroutine read_sfcelv_multi(sno,ix, iy, ndays, fsfcelv,nXX, nYY, wse)
! ===============================================
        implicit none
        integer,intent(IN)                               :: sno
        integer,dimension(sno),intent(IN)                :: ix,iy
        !integer,dimension(sno),intent(IN)                :: ix2,iy2
        integer,intent(IN)                               :: ndays
        integer,intent(IN)                               :: nXX, nYY
        character*128,intent(IN)                         :: fsfcelv
        real,dimension(sno,ndays),intent(OUT)            :: wse

        real,dimension(nXX,nYY)                          :: sfcelv
        integer                                          :: irec
        integer                                          :: isno, ix_sno,iy_sno !,,ix2_sno,iy2_sno


        open(12, file=fsfcelv, form='unformatted', access='direct', &
            & recl=4*nXX*nYY)
    
        do irec = 1,ndays
            read(12,rec=irec) sfcelv
            do isno = 1, sno
                ix_sno = ix(isno)+1
                iy_sno = iy(isno)+1
                ! iy1_sno = iy1(isno)+1
                ! iy2_sno = iy2(isno)+1
                if (ix_sno .GT. 0 .and. iy_sno .GT. 0 ) then
                wse(isno,irec)=sfcelv(ix_sno,iy_sno) !+ outflw(ix2_sno,iy2_sno)
                ! else
                !   wse(isno,irec)=outflw(ix1_sno,iy1_sno)
                end if 
                !print *, irec, wse(irec)
            end do
        enddo

        close(12)
        end subroutine read_sfcelv_multi
!***************************************************************************