# AMPT_small_systems

## How to get AMPT

This code is made to run on AMPT version "Ampt-v1.26t5-v2.26t5".  
You can get this version of AMPT from the following link:
http://myweb.ecu.edu/linz/ampt/

## Description of each file

### event_plane.C

This code is copied from rcf on Mar 2, 2017. Guaranteed could be run directly from rcf.  

Written in C++ to apply Event Plane method with three-sub events resolution. Three-sub events come from BBCS, FVTXS, FCTXN.  

Applied 6 ranges centrality cut.
  
Three-sub events method formula could be found in paper "Methods for analyzing anisotropic flow in relativistic nuclear collisions", equation 16. Or ask anyone in Nagle Lab.  
  
Function "processEvent" is used to compute Q vector and calculate resolution. The first part of numerator (of resolution formula) is filled in first bin of TProfile "res_comp". The second part of numerator is filled in second bin of TProfile "res_comp". The denominator is filled in the third bin.  
  
TProfile "v2s[i]" contains v2_pt histogram without modified by resolution. TProfile "res_comp[i]" contains three parts of resolution formula.  
  
It also contains N_ch distribution of BBCS (and FVTXS if need), eta distribution.  
  
To run this code, use the following command: root event_plane.C++