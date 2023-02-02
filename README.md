# HydroDA
  Data assimilation algorithm developed for global-scale data assimilation using a physically based emperical localization technique. Data assimilation method is local ensemble trasnformation Kalman filter (LETKF).
  
# HydroDA folder strucutre
    |-src
    |-out
    |-inp
    |-dat 
    |-img
    |-etc

   src : contains source codes

   out : output folder

   inp : folder containg input runoff forcing

   dat : folder containg river realted data

   img : folder containing codes for preapring images

   etc : folder containing code for miscellaneous 
   
********************
# Simple steps for HydroDA run
1. Download and compile CaMa-Flood (http://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/index.html) (Yamazaki et al 2011) model
2. Prepare observation data (e.g.: HydroWeb) : add the path to params.py file
3. Prepare physically-based empirical local patches (https://github.com/MenakaRevel/Empirical_LocalPatch) (Revel et al,. 2018, 2019) : add details to params.py file
4. Prepare runoff data in HydroDA/inp folder
5. Check river related data at HydroDA/dat folder {specially for anomaly or normalized assimilation}
6. Edit the settings at params.py and run_mool.sh at HydroDA/gosh
7. Compile HydroDA fortran codes using compile.sh
```bash
sh compile.sh "yes"
```
7. Edit the experiment name in run_mool.sh
8. Run run_mool.sh
```bash
qsub run_mool.sh
```
==========================================================
## Data assimilation codes (modules) mod_xxx.f90
data_assim.f90
 |- mod_letkf.f90
 |- mod_obser.f90
 |- mod_patch.f90
 |- mod_varxf.f90

1. mod_letkf.f90 - module for calculation analysis using LETKF (model independet code uses intel mkl library)
2. mod_obser.f90 - module for procedures on observation (reading observation, observation error etc.)
3. mod_patch.f90 - module for importing emperical local patches or other conventional local patches. (Emperical local patches: https://github.com/MenakaRevel/Empirical_LocalPatch.git)
4. mod_varxf.f90 - module for getting CaMa-Flood variables
+++++
x. data_assim.f90 - main program for data assimilation (uses mod_letkf, mod_obser, and mod_patch)
==========================================================
Reference:
1. Revel, M., Yamazaki, D., & Kanae, S. (2018). Model Based Observation Localization Weighting Function for Amazon Mainstream. Journal of Japan Society of Civil Engineers, Ser. B1 (Hydraulic Engineering), 74(5), I_157-I_162. https://doi.org/10.2208/jscejhe.74.5_I_157
2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
3. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 1–34. https://doi.org/10.1029/2020wr027876
4. Yamazaki, D., Kanae, S., Kim, H., & Oki, T. (2011). A physically based description of floodplain inundation dynamics in a global river routing model. Water Resources Research, 47(4), 1–21. https://doi.org/10.1029/2010WR009726
5. Revel, M., Zhou, X., Yamazaki, D., & Kanae, S. (2020).  Assimilation of transformed water surface elevation to improve river discharge estimation in a continental-scale river. Hydrology and Earth System Sciences, 27, 1-25, https://doi.org/10.5194/hess-27-1-2023, 2023