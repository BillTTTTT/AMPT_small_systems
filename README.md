# AMPT_small_systems

## How to get AMPT

This code is made to run on AMPT version "Ampt-v1.26t5-v2.26t5".  
You can get this version of AMPT from the following link:
http://myweb.ecu.edu/linz/ampt/

## Description of each file

### event_plane.C

There are two version of event_plane.C uploaded. 
One is "first commit". The other is "updated"

#### "First commit" version

This code is written in C++ to apply Event Plane method with three-sub events resolution. Three-sub events come from BBCS, FVTXS, FCTXN.  
  
Three-sub events method formula could be found in paper "Methods for analyzing anisotropic flow in relativistic nuclear collisions", equation 16. Or ask anyone in Nagle Lab.  
  
function "processEvent" is used to compute Q vector and calculate resolution. The first part of numerator (of resolution formula) is filled in first bin of TProfile "res_comp". The second part of numerator is filled in second bin of TProfile "res_comp". The denominator is filled in the third bin.  
  
Two TProfile will be written into output file. TProfile "v2s" contains v2_pt histogram without modified by resolution. TProfile "res_comp" contains three parts of resolution formula.  
  
#### "Update" version

This code is copied from rcf on Mar 2, 2017. Guaranteed could be run directly from rcf.  
  
This version contains everything in "First commit" version. But it applied centrality cut of 6 range.
  
It also contains N_ch distribution of BBCS (and FVTXS if need), eta distribution.