# HydroDA
  Data assimilation algorithem developed for global-scale data assimilation using a physically based emperical localization technique.
  
# HydroDA folder strucutre
    -src
    -out
    -inp
    -dat 
    -img

   src : contains source codes

   out : output folder
        -experiment_name
            -assim_out
            -cama_out
            -cama_in

   inp : input runoff forcing
        -data name
            -Roff

   dat : river realted data
        -river basin
        -HydroWeb station data

   img : codes for preapring images
   
********************
# Simple steps for HydroDA run
1. Download and compile CaMa-Flood (http://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/index.html) (Yamazaki et al 2011) model
2. Prepare observation data (e.g.: HydroWeb) : add the path to params file
3. Prepare runoff data in HydroDA/inp folder
4. Check river related data at HydroDA/dat folder {specially for anomaly or normalized assimilation}
5. Edit the settings at params.py and run_mool.sh at HydroDA/gosh
6. Compile HydroDA fortran codes using compile.sh
```bash
sh compile.sh "yes"
```
7. Edit the experiment name in run_mool.sh
8. Run run_mool.sh
```bash
qsub run_mool.sh
```
==========================================================

Reference:
1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 1–34. https://doi.org/10.1029/2020wr027876
2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829

