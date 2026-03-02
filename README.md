**Network Simulations**
1. 



**QC decisions**

* units with missing (or heightened) activity for certain blocks are assigned .5
* units with .2 < FR < .8 with decent metrics are retained, but can easily be dropped later
* example borderline rejections in eg_unit_rejections



**directory explanation**

results/ptID/

* raw
* records
* osort_mat
  * nsx2mat: raw data converted from .nsx to .mat
  * sorted_mats: .mat files after Osort
  * figs: figs generated after Osort
  * figs_clust: figs generated of only cluster (units) to make QC easier

note that pt 202509 and 20250717 dont have single units
