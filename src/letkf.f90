module letkf
!====================================================================================
! purpose: Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! References:
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
!
! Developed based on Takemasa Miyoshi
!====================================================================================
! created by Ikeshima & Menaka
! Menaka@IIS 2023
!====================================================================================
!$ use omp_lib
use common

implicit none

public

contains

subroutine letkf_core(ne,nobs,HEf,Rdiag,Rwgt,Yo,HXb,parm_infl,min_infl,infl_flg,T,errflg)
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     ne               : ensemble size                                          
!     nobs             : total number of observation assimilated at the point
!     HEf(nobs,ne)     : obs operator times fcst ens perturbations
!     Rdiag(nobs)      : observation error variance
!     Rwgt(nobs)       : localization weighting function
!     Yo(nobs)         : observations
!     HXb(nobs)        : simulations in observation space
!     parm_infl        : covariance inflation parameter
!     minfl            : minimum covariance inflation parameter
!     infl_flg         : infaltion flag (infl_flg = -1 : adaptive inflation, infl_flg > 1 : fixed inflation)      
!   OUTPUT
!     parm_infl        : updated covariance inflation parameter
!     T(ne,ne)         : transformation matrix
!     errflg           : error flag on making trasformation matrix
!     error flags:
!                2 = error in Eigon value decomposition
!                3 = error in calcualting square root
!=======================================================================
implicit none
integer,intent(in)           :: ne                      
integer,intent(in)           :: nobs
real(r_size),intent(in)      :: HEf(nobs,ne)
real(r_size),intent(in)      :: Rdiag(nobs)
real(r_size),intent(in)      :: Rwgt(nobs)
real(r_size),intent(in)      :: Yo(nobs)
real(r_size),intent(in)      :: HXb(nobs)
real(r_size),intent(inout)   :: parm_infl
real(r_size),intent(in)      :: min_infl
real(r_size),intent(in)      :: infl_flg
real(r_size),intent(out)     :: T(ne,ne)
integer(r_size),intent(out)  :: errflg

real(r_size)                 :: Rinv(nobs,nobs),dep(nobs)
real(r_size)                 :: UNI(ne,ne),HETRHE(ne,ne),VDVT(ne,ne),la_p(ne),U_p(ne,ne),la(ne),U(ne,ne)
real(r_size)                 :: Dinv(ne,ne),Dsqr(ne,ne),Pa(ne,ne),Pasqr(ne,ne)!,Tvec(ne,ne)
real(r_size)                 :: Tvec(ne)
real(r_size),allocatable     :: work(:),iwork(:),ifail(:),isuppz(:)
integer                      :: i,j,m
real(r_size)                 :: rho,rho_fixed !,min_infl ! covariance inflation parameters
integer                      :: lwork,liwork ! added for large ensemble size
integer                      :: info,info2
real(r_size)                 :: HEfR(ne,nobs),HPH(nobs,nobs) !HEf(nobs,ne),
real(r_size)                 :: parm(4)
real                         :: sigma_o,gain,sigma_b ! background variance
! real(r_size),allocatable     :: U(:,:),la(:,:),Dinv(:,:),Dsqr(:,:)

! real(r_size)                 :: W(ne,ne),Pa(ne,ne),Pasqr(ne,ne)!,HETRHE(ne,ne)!,HPH(:,:)
! real(r_size)                 :: rho,rho_fixed,min_infl ! covariance inflation parameters
! real(r_size)                 :: VDVTmax!,xf_m55,xa_m55,xt_m55,xf_err,xa_err,delta,traceR,traceHPH
! real(r_size)                 :: la(ne),U(ne,ne),Dinv(ne,ne),VDVT(ne,ne),Dsqr(ne,ne) !,yo(:),Ef(:,:),la_p(:),U_p(:,:)


! initialize error flag
errflg=0

! observation departure (Yo-HXb)
dep=Yo-HXb

! R inverse
Rinv=0.
do i=1,nobs
    Rinv(i,i)=(Rdiag(i)/(Rwgt(i)+1e-20))**(-1.)
end do


!===========================================
! covariance inflation parameter
! rho=1.0d0
if (infl_flg == -1.0) then
    rho=parm_infl
    ! if (parm_infl > 1.0d0) then
    if (parm_infl<min_infl) rho=min_infl
else
    rho=infl_flg
end if

!===========================================
! calculate VDVT
UNI=0
UNI=RESHAPE([(1,(0,i=1,ne),j=1,ne-1),1],[ne,ne])
HETRHE=matmul(matmul(TRANSPOSE(HEf),Rinv),HEf)
! VDVTmax=maxval(abs(HETRHE))
VDVT=real(ne-1.)*UNI/rho+HETRHE

!===========================================
! Eigon value decomposition
! Diagonalized VDVT to calculate inverse matrix
la_p=0
U_p=0
lwork=max(1,26*ne)
liwork=max(1,10*ne)
allocate(work(lwork),iwork(lwork),ifail(lwork),isuppz(lwork))
!call ssyevx("V","A","U",ens_num,VDVT,ens_num,1e-5,1e5,1,2,1.2e-38*2.,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)
!call ssyevx("V","A","U",ens_num,VDVT,ens_num,-1e20,1e20,1,ens_num,-1.0,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)
!call ssyevx("V","I","U",ens_num,VDVT,ens_num,-1e20,1e20,1,ens_num,-1.0,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)
!call ssyevr("V","A","U",ens_num,VDVT,ens_num,-1e-20,1e20,1,ens_num,2.0*2.3e-38,m,la_p,U_p,ens_num,isuppz,work,1000,iwork,1000,info)
! call ssyevr("V","A","U",ens_num,VDVT,ens_num,-1e20,1e20,1,ens_num,-1.0,m,la_p,U_p,ens_num,isuppz,work,1000,iwork,1000,info)
call ssyevr("V","A","U",ne,VDVT,ne,-1e20,1e20,1,ne,-1.0,m,la_p,U_p,ne,isuppz,work,lwork,iwork,lwork,info)

!===========================================
! calculate transforamtion matrix
if (m<ne) then
    errflg=2
    T=0.0
else
    ! allocate(U(m,m),la(m))
    U=0
    la=0
    do i=1,ne
        do j=1,ne
            U(i,j)=U_p(i,j)
        end do
        la(i)=la_p(i)
    end do
    ! allocate(Dinv(m,m),Dsqr(m,m))
    Dinv=0
    Dsqr=0
    ! calc Dinv,Dsqr
    if(info==0)then
        do i=1,ne
            Dinv(i,i)=(la(i)+1e-20)**(-1)
        end do
        !===========================================
        ! calculate square root
        !=write(78,*) "Dinv",Dinv
        Dsqr=Dinv
        info2=0
        call spotrf("U",m,Dsqr,m,info2)
        if(info2/=0) then
            errflg=3
            T=0.0
        end if
    end if
    ! Pa    = VD-1VT
    ! Pasqr = VD-1/2VT
    Pa=0
    Pasqr=0
    Tvec=0
    Pa   =matmul(matmul(U,Dinv),transpose(U))
    Pasqr=matmul(matmul(U,Dsqr),transpose(U))
    Tvec = matmul(matmul(matmul(Pa,transpose(HEf)),Rinv),dep)
    do i = 1,ne
        T(:,i) = Tvec(i) + sqrt(ne-1.)*Pasqr(:,i)
    end do
end if

!===========================================
!-- adapative inflation parameter estimation
parm=0.0d0
! allocate(dep(ovs),HEf(ovs,ens_num),HEfR(ens_num,ovs),HPH(ovs,ovs))
!write(*,*) yo-matmul(H,xf_m)
! dep=yo-matmul(H,xf_m)
! HEf=matmul(H,Ef)
HEfR=matmul(TRANSPOSE(HEf),Rinv)
HPH=matmul(HEf,TRANSPOSE(HEf))
do i=1,nobs
    parm(1)=parm(1)+dep(i)*dep(i)*Rdiag(i)**(-1.)
end do
!do num=1,ens_num
do i=1,nobs
    parm(2)=parm(2)+HPH(i,i)*Rdiag(i)**(-1.)
end do
!end do
parm(2)=parm(2)/real(ne-1,r_size)
parm(3)=sum(Rwgt)
parm(4)=(parm(1)-parm(3))/(parm(2)+1.0e-20) - rho
sigma_o=2.0d0/(parm(3)+1.0e-20) * ((rho*parm(2) + parm(3))/(parm(2)+1.0e-20))**2
gain=sigma_b**2 / (sigma_o + sigma_b**2)
! write(73,*) "+++++++++++++++++++++++++++++++++++"
! write(73,*) gain, parm(4),parm(1),parm(2),parm(3)
rho=rho+ gain* parm(4)
! write(*,*)"rho",rho,rho_min
if (rho<min_infl) rho=min_infl
parm_infl=rho
! write(73,*) lon_cent,lat_cent,"rho",rho,"rho_min",rho_min

return
deallocate(work,iwork,ifail,isuppz)

end subroutine letkf_core

end module letkf