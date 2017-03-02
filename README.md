# AMPT_small_systems

## How to get AMPT

This code is made to run on AMPT version "Ampt-v1.26t5-v2.26t5".  
You can get this version of AMPT from the following link:
http://myweb.ecu.edu/linz/ampt/

## Description of each file

### event_plane.C

This code is directly copied from rcf on Mar 2, 2017. Guaranteed could be run from rcf.  

Written in C++, applying Event Plane method with three-sub events resolution. Three-sub events come from BBCS, FVTXS, FCTXN.  

Applied 6 ranges centrality cut.
  
Three-sub events method formula could be found in paper "Methods for analyzing anisotropic flow in relativistic nuclear collisions", equation 16. Or ask anyone in Nagle Lab.  
  
Function "processEvent" is used to compute Q vector and calculate resolution. The first part of numerator (of resolution formula) is filled in first bin of TProfile "res_comp". The second part of numerator is filled in second bin of TProfile "res_comp". The denominator is filled in the third bin.  
  
TProfile "v2s[i]" contains v2_pt histogram without modified by resolution. TProfile "res_comp[i]" contains three parts of resolution formula.  
  
It also contains N_ch distribution of BBCS (and FVTXS if need), eta distribution.  
  
To run this code, use the following command: root event_plane.C++  

### cumulant.C

This code is directly copied from rcf on Mar 2, 2017. Guaranteed could be run from rcf.

Written in C++, applying Multi-particle Correlation (Cumulant) Method with Ron's help with 4 and 6 particle cumulant equation. This file contains reference flow of 2,4,6 particle cumulant (usually used part in my project). And diferential flow of 2,4 particle cumulant (not often used in my projects). Formula and theory could be found from following paper: 

Flow analysis with cumulants: direct calculations, A. Bilandzic, R. Snellings, S. Voloshin, Phys.Rev.C83:044913,2011, (or arXiv:1010.0233v2)  
  
Also applied an eta gap in 2 particle cumulant calculation.  
  
#### Functions 
  
"def_ave_2particle_with_gap" 			to calculate <2'> term in diferential flow, with eta gap.  

"def_ave_2particle_correlation"		to calculate <2'> term in diferential flow, without eta gap.  
  
"def_ave_4particle_correlation"		to calculate <4'> term in diferential flow.  
  
"with_gap_calculation"				to calculate <2> term in reference flow, with eta gap.  
  
"ave_2particle_correlation" 			to calculate <2> term in reference flow, without eta gap.  
  
"ave_4particle_correlation" 			to calculate <4> term in reference flow.  

"ave_6particle_correlation" 			to calculate <6> term in reference flow.  

The variable name "raai" (i.e. raa2) means reference flow: <<i>> (i.e. <<2>>). "daai" means diferential flow <<i>>.  
  
#### Output
  
"comp_Ncharge"			Multiplicity dependence of 2 particle cumulant method with eta gap.  

"raa2_Ncharge"			Multiplicity dependence of <<2>> in reference flow.  
  
"raa4_Ncharge"			Multiplicity dependence of <<4>> in reference flow.  
  
"raa6_Ncharge"			Multiplicity dependence of <<6>> in reference flow.  
  
"daa2_Ncharge"			Multiplicity dependence of <<2'>> in differential flow.  
  
"daa4_Ncharge"			Multiplicity dependence of <<4'>> in differential flow.  
  
"dnch" 					N_ch distribution





















