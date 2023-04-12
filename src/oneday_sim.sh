#!/bin/sh
#==========================================================
# CaMa-Flood sample go script (1) global 15min simulation
# -- Multi 1-year simulations (2000 spinup -> 2000 -> 2001)
# -- Daily runoff forcing (plain binary) at 1deg resolution
#
# (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is 
#  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and limitations under the License.
#==========================================================

#*** PBS setting when needed
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -V
#========
#edited by Menaka@IIS for SWOTDA 2019/11/05
#================================================
# input settings
orgDIR=`pwd`

in_year=$1
in_month=$2
in_date=$3

ens_num=$4

CAMADIR=$5

looptype=$6

cpunums=$7

runname=$8

EXP_DIR=$9

mapname=${10}

cal=${11}

DA_dir=${12}

# echo "looptype" $looptype
#=================================================
#cd ${CAMADIR}/gosh
#ens_num=$(printf '%03d' $(($runens*$manens)))
# pwd
# echo ${CAMADIR}
echo $EXP_DIR
# years,months,dates in arranged digit
ar_year=$in_year
ar_month=$in_month
ar_date=$in_date
#----
if [ $looptype = "true" ];then
	ar_ens="T"$ens_num
elif [ $looptype = "open" ];then
	ar_ens="C"$ens_num
else
	ar_ens="A"$ens_num
fi
#================================================
# (0) Basic Setting (for workstation)

#*** 0a. Set CaMa-Flood base directory
BASE=$CAMADIR                         #   CaMa-Flood directory
OUTBASE="$EXP_DIR/CaMa_out"    		   #   base for output => CaMa_out
INBASE="../../CaMa_in"                #   base for input => CaMa_in

# echo $BASE

#*** 0b. Set dynamic library if needed
export IFORTLIB="/opt/intel/lib:/opt/intel/mkl/lib"
export DYLD_LIBRARY_PATH="${IFORTLIB}:${DYLD_LIBRARY_PATH}"

#*** 0c. OpenMP thread number
export OMP_NUM_THREADS=$cpunums       # OpenMP cpu num

#================================================
# (1) Experiment setting
# -- some non-default options can be modified in NAMELIST section

#============================
#*** 1a. Experiment directory setting
EXP=$ar_year$ar_month$ar_date$ar_ens        # experiment name (output directory name)
RDIR=${OUTBASE}/${EXP}                      # directory to run CaMa-Flood
EXE="MAIN_cmf"                              # Execute file name
PROG=${BASE}/src/${EXE}                     # location of Fortran main program
NMLIST="./input_cmf.nam"                    # standard namelist
LOGOUT="./log_CaMa.txt"                     # standard log output


#============================
#*** 1b. Model physics option
DT=86400                                    # base DT (modified in physics loop by LADPSTP)
LADPSTP=".TRUE."                            # .TRUE. for adaptive time step

LFPLAIN=".TRUE."                            # .TRUE. to activate floodplain storage
LKINE=".FALSE."                             # .TRUE. to use kinematic wave equation
LFLDOUT=".TRUE."                            # .TRUE. to activate floodplain discharge
LPTHOUT=".TRUE."                            # .TRUE. to activate bifurcation flow, mainly for delta simulation
LDAMOUT=".FALSE."                           # .TRUE. to activate reservoir operation (under development)


#============================
#*** 1c. simulation time
# simulation for spinup is 1 day
#YSTA=2004                                   # start year ( from YSTA / Jan  1st _ 00:00)
#YEND=2004                                   # end   year (until YEND / Dec 31st _ 24:00)
SPINUP=1                                     # [0]: zero-storage start, [1]: from restart file
NSP=0                                        # spinup repeat time


#============================
#*** 1d. spinup setting

#* input restart file
LRESTART=".TRUE." # see (3) set each year   # TRUE. to use restart initial condition
#CRESTSTO="" # see (3) set each year         # input restart FIle
if [ $looptype = "true" ];then
	CRESTSTO=$INBASE"/restart/true/restart"$ar_year$ar_month$ar_date"T000.bin" #restart file name
elif [ $looptype = "open" ];then
	CRESTSTO=$INBASE"/restart/open/restart"$ar_year$ar_month$ar_date"C"$ens_num".bin" #restart file name
else
	CRESTSTO=$INBASE"/restart/assim/restart"$ar_year$ar_month$ar_date"A"$ens_num".bin" #restart file name
fi
# echo $CRESTSTO
LSTOONLY=".TRUE."                          # .TRUE. for storage only restart (for assimilation)


#* output restart file
CRESTDIR="./"                               # output restart file directory
CVNREST="restart"                           # output restart file prefix
LRESTCDF=".FALSE."                          # .TRUE. to use netCDF restart file
IFRQ_RST="0"                                # output restat frequency.
                                            # [0]: only at last time, [1,2,3,...,24] hourly restart, [30]: monthly restart


#============================
#*** 1e. forcing setting
IFRQ_INP="24"                               # input forcing frequency: [1,2,3,...,24] hour
DROFUNIT="86400000"        # [mm/day->m/s]  # runoff unit conversion
if [ $runname = "E2O" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ECMWF000" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ECMWF050" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA20CM" ];then
     DROFUNIT="1000"       # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ELSE_KIM2009" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "VIC_BC" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "isimip3a" ];then
     DROFUNIT="1000"       # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA5" ];then
     DROFUNIT="86400"      # [m/day->m/s]   # runoff unit conversion
fi

#----- for plain binary runoff forcing
LINPCDF=".FALSE."                           # true for netCDF runoff
LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
#CROFDIR="${BASE}/inp/test_1deg/runoff/"     # runoff directory
#CROFPRE="Roff____"                          # runoff prefix/suffix  
#CROFSUF=".one"                              #   $(CROFPRE)YYYYMMDD$(CROFSUF)
#CROFDIR="/cluster/data6/menaka/SWOTDA_vRunoff_EmpPatch/CaMa_in/E2O/Roff"     # runoff directory
CROFPRE="Roff__"                          # runoff prefix/suffix
#CROFSUF="001.one"
if [ $looptype = "true" ];then
	CROFDIR="${INBASE}/CaMa_in/$runname/Roff_TRUE"    #   runoff directory
	CROFSUF="T000.one" 
else
	CROFDIR="${INBASE}/$runname/Roff"    #   runoff directory
	CROFSUF="${ens_num}.one" 
fi

###** sub-surface runoff scheme (not available with plain binary runoff)
LROSPLIT=".FALSE."                          # .TRUE. for sub-surface runoff
###CSUBDIR="NONE"                              # sub-surface runoff directory
###CSUBPRE="NONE"                              # sub-surface runoff prefix/suffix  
###CSUBSUF="NONE"                              #   $(PREFIX)YYYYMMDD$(SUFFIX)

#----- for netCDF runoff forcing ###
###LINPCDF=".TRUE."                              # true for netCDF runoff
###LINTERP=".TRUE."                              # .TRUE. to interporlate with input matrix
###LINTERPCDF=".FALSE."                          # .TRUE. to use netCDF input matrix
###CROFDIR="${BASE}/inp/test_15min_nc/"          # runoff directory
###CROFPRE="e2o_ecmwf_wrr2_glob15_day_Runoff_"   # runoff prefix/suffix  
###CROFCDF=""     # see (3) set each year        # netCDF runoff file
###CVNROF="Runoff"                               # netCDF runoff    variable name
###CVNSUB=""                                     # netCDF runoffsub variable name
###SYEARIN=""     # see (3) set each year        #   netCDF runoff file, start date
###SMONIN=""      # see (3) set each year
###SDAYIN=""      # see (3) set each year
###SHOURIN=""     # see (3) set each year


#============================
#*** 1f. river map & topography
FMAP="${BASE}/map/${mapname}"                # map directory
CDIMINFO="${FMAP}/diminfo_test-1deg.txt"    # dimention information file
CINPMAT="${FMAP}/inpmat_test-1deg.bin"        # runoff input matrix for interporlation
#CDIMINFO="${FMAP}/diminfo_test-15min_nc.txt" # dimention information file
#CINPMAT=${FMAP}/inpmat_test-15min_nc.bin     # runoff input matrix for interporlation
#CDIMINFO="${FMAP}/diminfo_test-15min.txt" # dimention information file
#CINPMAT=${FMAP}/inpmat_test-15min.bin     # runoff input matrix for interporlation
if [ $runname = "E2O" ] ; then
     CDIMINFO="${FMAP}/diminfo-15min.txt" # dimention information file
     CINPMAT="${FMAP}/inpmat-15min.bin"     # runoff input matrix for interporlation
elif [ $runname = "ECMWF000" ];then
    CDIMINFO="${FMAP}/diminfo-15min.txt" # dimention information file
    CINPMAT="${FMAP}/inpmat-15min.bin"   # runoff input matrix for interporlation
elif [ $runname = "ECMWF050" ] ; then
     CDIMINFO="${FMAP}/diminfo-15min.txt" # dimention information file
     CINPMAT="${FMAP}/inpmat-15min.bin"     # runoff input matrix for interporlation
elif [ $runname = "ERA20CM" ] ; then
     CDIMINFO="${FMAP}/diminfo-1deg.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-1deg.bin"      # runoff input matrix for interporlation
elif [ $runname = "ELSE_KIM2009" ] ; then
     CDIMINFO="${FMAP}/diminfo-1deg.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-1deg.bin"      # runoff input matrix for interporlation
elif [ $runname = "VIC_BC" ] ; then
     CDIMINFO="${FMAP}/diminfo-15min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-15min.bin"      # runoff input matrix for interporlation
elif [ $runname = "isimip3a" ] ; then
     CDIMINFO="${FMAP}/diminfo-30min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-30min.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "ERA5" ] ; then
     CDIMINFO="${FMAP}/diminfo-06min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-06min.bin"      # runoff input matrix for interporlation
fi

#----- for plain binary map input
#** basic topography
LMAPCDF=".FALSE."                           # .TRUE. for netCDF map
CNEXTXY="${FMAP}/nextxy.bin"                # downstream xy (river network map)
CGRAREA="${FMAP}/ctmare.bin"                # unit-catchment area   [m2]
CELEVTN="${FMAP}/elevtn.bin"                # channel top elevation [m]
CNXTDST="${FMAP}/nxtdst.bin"                # downstream distance   [m]
CRIVLEN="${FMAP}/rivlen.bin"                # channel length        [m]
CFLDHGT="${FMAP}/fldhgt.bin"                # floodplain elevation profile (height above 'elevtn') [m]

#** channel parameter
###CRIVWTH=${FMAP}/rivwth.bin"              # channel width [m] (empirical power-low)
CRIVWTH="${FMAP}/rivwth_gwdlr.bin"          # channel width [m] (GWD-LR + filled with empirical)
CRIVHGT="${FMAP}/rivhgt.bin"                # channel depth [m] (empirical power-low)
if [ $cal = "yes" ];then
  CRIVHGT="${FMAP}/rivhgt_Xudong.bin"         # channel depth [m] (Xudong et al 2021)
elif [ $cal = "corrupt" ];then
  CRIVHGT="${FMAP}/rivhgt_corrupt.bin"         # channel depth [m] (Corrupted rivhgt)
fi
CRIVMAN="${FMAP}/rivman.bin"                # manning coefficient river (The one in flood plain is a global parameter; set $PMANFLD below.)
#if [ $looptype = "true" ] ; then
#    CRIVMAN="${INBASE}/assim_out/rivman/rivmanTRUE.bin"
#    #CRIVMAN="${FMAP}/rivmanTRUE.bin"
#else
#    CRIVMAN="${INBASE}/assim_out/rivman/rivmanCORR.bin"
#    #CRIVMAN="${FMAP}/rivmanCORR.bin"
#fi
# echo $CRIVMAN

#** bifurcation channel info
CPTHOUT="${FMAP}/bifprm.txt"                #   bifurcation channel list

###** groundwater delay (not available in plain binary runoff/map)
LGDWDLY=".FALSE."                           # .TRUE. to actuvate groundwater delay
#CGDWDLY=""                                 # ground water delay map

###** mean sea level
LMEANSL=".FALSE."                           # .TRUE. to use mean sea level data
#CMEANSL=""                                 # mean sea level map

#----- for netCDF map input 
###LMAPCDF=".TRUE."                         # .TRUE. for netCDF map
###CRIVCLINC=""                             # netCDF topography map
###CRIVPARNC=""                             # netCDF river parameters
###CMEANSLNC=""                             # netCDF mean sea level


#============================
#*** 1g. Dynamic Boundary Sea Level (not default)
LSEALEV=".FALSE."                           # .TRUE. to activate dynamic sea level 
###LSEALEVCDF=".FALSE."
###CSEALEVDIR="NONE"                        # Sea level boundary DIRECTORY
###CSEALEVPRE="NONE"                        # Sea level boundary PREFIX
###CSEALEVSUF="NONE"                        # Sea level boundary SUFFIX
###CSEALEVCDF="NONE"                        # * Sea level netCDF file name
###CVNSEALEV="sealev"                       # * Sea Level netCDF variable name
###SYEARSL=1                                # * netCDF sea level start year
###SMONSL=1                                 # * netCDF sea level start year
###SDAYSL=1                                 # * netCDF sea level start year
###SHOURSL=0                                # * netCDF sea level start year
###NSTATIONS=1                              # sea level data points
###CSLMAP="NONE"                            # sea level sta->XY conversion table

#============================
#*** 1h. Output Settings 
LOUTPUT=".TRUE."                            # .TRUE. to use CaMa-Flood standard output
IFRQ_OUT=24                                 # output frequency: [1,2,3,...,24] hour

LOUTCDF=".FALSE."                           # .TRUE. netCDF output, .FALSE. plain binary output
COUTDIR="./"                                # output directory 
#CVARSOUT="outflw,storge,fldfrc,maxdph,flddph" # list output variable (comma separated)
#CVARSOUT="rivout,rivsto,rivdph,rivvel,fldout,fldsto,flddph,fldfrc,fldare,sfcelv,outflw,storge,pthflw,pthout,maxsto,maxflw,maxdph" # list output variable (comma separated)
#CVARSOUT="rivout,rivsto,rivdph,fldout,flddph,fldfrc,fldare,sfcelv,outflw,storge,maxsto,maxflw,maxdph" # list output variable (comma separated)
CVARSOUT="rivout,rivsto,fldout,flddph,fldfrc,fldare,sfcelv,outflw,storge,maxsto,maxflw,maxdph" # list output variable (comma separated)
COUTTAG=""  # see (3) set each year         #   output tag $(COUTDIR)/$(VARNAME)$(OUTTAG).bin

##### Model Parameters ################
PMANRIV="0.03D0"                            # manning coefficient river
PMANFLD="0.10D0"                            # manning coefficient floodplain
PCADP="0.7"                                 # satety coefficient for CFL condition
PDSTMTH="10000.D0"                          # downstream distance at river mouth [m]



#================================================
# (2) Initial setting

#*** 2a. create running dir 
mkdir -p ${RDIR}
cd ${RDIR}

#*** 2b. for new simulation, remove old files in running directory

if [ ${SPINUP} -eq 0 ]; then
  rm -rf ${RDIR}/????-sp*
  rm -rf ${RDIR}/*.bin
  rm -rf ${RDIR}/*.pth
  rm -rf ${RDIR}/*.vec
  rm -rf ${RDIR}/*.nc
  rm -rf ${RDIR}/*.log
  rm -rf ${RDIR}/*.txt
  rm -rf ${RDIR}/restart*
else
  NSP=0  # restart, no spinup
fi

#================================================
# (3) For each simulation year, modify setting
#--  loop 1-day simulation from $YSTART to $YEND
YSTA=$in_year
ISP=1           ## spinup count
IYR=${YSTA}     ## curent year
CYR=`printf %04d ${IYR}`

#*** 3a. modify restart setting
#if [ ${SPINUP} -eq 0 ];then
#  LRESTART=".FALSE."                  ## from zero storage
#  CRESTSTO=""
#else
#  LRESTART=".TRUE."
#  CRESTSTO=$CRESTSTO
##    CRESTSTO="${CVNREST}${CYR}010100.bin"   ## from restart file
##    CRESTSTO="${CVNREST}${CYR}010100.nc"    ## from restart file
#fi

#*** 3b. update start-end year
SYEAR=$in_year
SMON=$in_month
SDAY=$in_date
SHOUR=0

EYEAR=`python ${DA_dir}/src/calc_end_date.py $in_year $in_month $in_date "year"`
EMON=`python ${DA_dir}/src/calc_end_date.py $in_year $in_month $in_date "month"`
EDAY=`python ${DA_dir}/src/calc_end_date.py $in_year $in_month $in_date "date"`
EHOUR=0

ln -sf $PROG $EXE

#*** 3c. update input / output file data
CSYEAR=`printf %04d ${SYEAR}`
COUTTAG=${CSYEAR}                  # output file tag

#CROFCDF="${CROFDIR}/${CROFPRE}${CSYEAR}.nc"  # input netCDF runoff file
#SYEARIN=$IYR
#SMONIN=1
#SDAYIN=1
#SHOURIN=0

#================================================
# (4) Create NAMELIST for simulation year
# it is OK to remove optional variables (set to default in CaMa-Flood)

rm -f ${NMLIST}

#*** 0. config
cat >> ${NMLIST} << EOF
&NRUNVER
LADPSTP  = ${LADPSTP}                  ! true: use adaptive time step
LFPLAIN  = ${LFPLAIN}                  ! true: consider floodplain (false: only river channel)
LKINE    = ${LKINE}                    ! true: use kinematic wave
LFLDOUT  = ${LFLDOUT}                  ! true: floodplain flow (high-water channel flow) active
LPTHOUT  = ${LPTHOUT}                  ! true: activate bifurcation scheme
LDAMOUT  = ${LDAMOUT}                  ! true: activate dam operation (under development)
LROSPLIT = ${LROSPLIT}                 ! true: input if surface (Qs) and sub-surface (Qsb) runoff
LGDWDLY  = ${LGDWDLY}                  ! true: Activate ground water reservoir and delay
LSLPMIX  = .FALSE.                     ! true: activate mixed kinematic and local inertia based on slope
LMEANSL  = ${LMEANSL}                  ! true: boundary condition for mean sea level
LSEALEV  = ${LSEALEV}                  ! true: boundary condition for variable sea level
LRESTART = ${LRESTART}                 ! true: initial condition from restart file
LSTOONLY = ${LSTOONLY}                 ! true: storage only restart (mainly for data assimilation)
LOUTPUT  = ${LOUTPUT}                  ! true: use standard output (to file)
LGRIDMAP = .TRUE.                      ! true: for standard XY gridded 2D map
LLEAPYR  = .TRUE.                      ! true: neglect leap year (Feb29 skipped)
LMAPEND  = .FALSE.                     ! true: for map data endian conversion
LBITSAFE = .FALSE.                     ! true: for Bit Identical simulation (avoid OSM ATOMIC)
/
&NDIMTIME
CDIMINFO = "${CDIMINFO}"               ! text file for dimention information
DT       = ${DT}                       ! time step length (sec)
IFRQ_INP = ${IFRQ_INP}                 ! input forcing update frequency (hour)
/
&NPARAM
PMANRIV  = ${PMANRIV}                  ! manning coefficient river
PMANFLD  = ${PMANFLD}                  ! manning coefficient floodplain
PGRV     = 9.8D0                       ! gravity accerelation
PDSTMTH  = ${PDSTMTH}                  ! downstream distance at river mouth [m]
PCADP    = ${PCADP}                    ! CFL coefficient
PMINSLP  = 1.D-5                       ! minimum slope (kinematic wave)
IMIS     = -9999                       ! missing value for integer
RMIS     = 1.E20                       ! missing value for real*4
DMIS     = 1.E20                       ! missing value for real*8
CSUFBIN  = '.bin'                      ! file suffix for plain binary 2D map
CSUFVEC  = '.vec'                      ! file suffix for plain binary 1D vector
CSUFPTH  = '.pth'                      ! file suffix for plain binary bifurcation channel
CSUFCDF  = '.nc'                       ! file suffix for netCDF
/
EOF

#*** 1. time
cat >> ${NMLIST} << EOF
&NSIMTIME
SYEAR   = ${SYEAR}                     ! start year
SMON    = ${SMON}                      !  month
SDAY    = ${SDAY}                      !  day
SHOUR   = ${SHOUR}                     !  hour
EYEAR   = ${EYEAR}                     ! end year
EMON    = ${EMON}                      !  month
EDAY    = ${EDAY}                      !  day
EHOUR   = ${EHOUR}                     !  hour
/
EOF

#*** 2. map
cat >> ${NMLIST} << EOF
&NMAP
LMAPCDF    = ${LMAPCDF}                ! * true for netCDF map input
CNEXTXY    = "${CNEXTXY}"              ! river network nextxy
CGRAREA    = "${CGRAREA}"              ! catchment area
CELEVTN    = "${CELEVTN}"              ! bank top elevation
CNXTDST    = "${CNXTDST}"              ! distance to next outlet
CRIVLEN    = "${CRIVLEN}"              ! river channel length
CFLDHGT    = "${CFLDHGT}"              ! floodplain elevation profile
CRIVWTH    = "${CRIVWTH}"              ! channel width
CRIVHGT    = "${CRIVHGT}"              ! channel depth
CRIVMAN    = "${CRIVMAN}"              ! river manning coefficient
CPTHOUT    = "${CPTHOUT}"              ! bifurcation channel table
CGDWDLY    = "${CGDWDLY}"              ! Groundwater Delay Parameter
CMEANSL    = "${CMEANSL}"              ! mean sea level
CRIVCLINC  = "${CRIVCLINC}"            ! * river map netcdf
CRIVPARNC  = "${CRIVPARNC}"            ! * river parameter netcdf (width, height, manning, ground water delay)
CMEANSLNC  = "${CMEANSLNC}"            ! * mean sea level netCDF
/
EOF

#*** 3. restart
cat >> ${NMLIST} << EOF
&NRESTART
CRESTSTO = "${CRESTSTO}"               ! restart file
CRESTDIR = "${CRESTDIR}"               ! restart directory
CVNREST  = "${CVNREST}"                ! restart variable name
LRESTCDF = ${LRESTCDF}                 ! * true for netCDF restart file
IFRQ_RST = ${IFRQ_RST}                 ! restart write frequency (1-24: hour, 0:end of run)
/
EOF

#*** 4. forcing
if [ ${LINPCDF} = ".FALSE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
LINPEND  = .FALSE.                     ! true for runoff endian conversion
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFDIR  = "${CROFDIR}"                ! runoff             input directory
CROFPRE  = "${CROFPRE}"                ! runoff             input prefix
CROFSUF  = "${CROFSUF}"                ! runoff             input suffix
CSUBDIR  = "${CSUBDIR}"                ! sub-surface runoff input directory
CSUBPRE  = "${CSUBPRE}"                ! sub-surface runoff input prefix
CSUBSUF  = "${CSUBSUF}"                ! sub-surface runoff input suffix
/
EOF

elif [ ${LINPCDF} = ".TRUE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
LINPEND  = .FALSE.                     ! true for runoff endian conversion
LITRPCDF = ${LINTERPCDF}               ! * true for netCDF input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFCDF  = "${CROFCDF}"                ! * netCDF input runoff file name
CVNROF   = "${CVNROF}"                 ! * netCDF input runoff variable name
CVNSUB   = "${CVNSUB}"                 ! * netCDF input runoff variable name
SYEARIN  = ${SYEARIN}                  ! * netCDF input start year
SMONIN   = ${SMONIN}                   ! * netCDF input start year
SDAYIN   = ${SDAYIN}                   ! * netCDF input start year
SHOURIN  = ${SHOURIN}                  ! * netCDF input start year
/
EOF

fi # (if LINPCDF)

#*** 5. outputs
cat >> ${NMLIST} << EOF
&NOUTPUT
COUTDIR  = "${COUTDIR}"                ! OUTPUT DIRECTORY
CVARSOUT = "${CVARSOUT}"               ! Comma-separated list of output variables to save 
COUTTAG  = "${COUTTAG}"                ! Output Tag Name for each experiment
LOUTVEC  = .FALSE                      ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
LOUTCDF  = ${LOUTCDF}                  ! * true for netcdf outptu false for binary
NDLEVEL  = 0                           ! * NETCDF DEFLATION LEVEL 
IFRQ_OUT = ${IFRQ_OUT}                 ! output data write frequency (hour)
/
EOF

#### 6. sea level (optional) 
#cat >> ${NMLIST} << EOF
#&NBOUND
#LSEALEVCDF =  ${LSEALEVCDF}            ! * true : netCDF sea level boundary
#CSEALEVDIR = "${CSEALEVDIR}"           ! Sea level boundary DIRECTORY
#CSEALEVPRE = "${CSEALEVPRE}"           ! Sea level boundary PREFIX
#CSEALEVSUF = "${CSEALEVSUF}"           ! Sea level boundary SUFFIX
#CSEALEVCDF = "${CSEALEVCDF}"           ! * Sea level netCDF file name
#CVNSEALEV  = "${CVNSEALEV}"            ! * Sea Level netCDF variable name
#SYEARSL    = ${SYEARSL}                ! * netCDF sea level start year
#SMONSL     = ${SMONSL}                 ! * netCDF sea level start year
#SDAYSL     = ${SDAYSL}                 ! * netCDF sea level start year
#SHOURSL    = ${SHOURSL}                ! * netCDF sea level start year
#NSTATIONS  = ${NSTATIONS}              ! sea level data points
#CSLMAP     = "${CSLMAP}                ! station to XY conversion table
#IFRQ_SL    = ${IFRQ_SL}                ! sea level boundary update frequency (min)
#/
#EOF

#================================================
# (5) Execute main program

echo "start: ${SYEAR}" `date`  >> log.txt
time ./${EXE}                  >> log.txt 
echo "end:   ${SYEAR}" `date`  >> log.txt

mv ${LOGOUT} log_CaMa-${CYR}.txt

# move restart file
RYR=`printf %04d ${EYEAR}`
RMO=`printf %02d ${EMON}`
RDA=`printf %02d ${EDAY}`
cp restart${RYR}${RMO}${RDA}00.bin restart${RYR}${RMO}${RDA}.bin

###################
# going back to original directory
cd $orgDIR

exit 0
