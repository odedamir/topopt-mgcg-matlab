# topopt-mgcg-matlab
MATLAB codes for efficient 3-D topology optimization
Codes use multigrid CG, approximate sensitivity analysis, recycling precoditioners

This repository includes matlab codes that were used in the following papers:
1) Amir, O., Aage, N., & Lazarov, B. S. (2014). On multigrid-CG for efficient topology optimization. 
Structural and Multidisciplinary Optimization, 49(5), 815-829.
2) Amir, O. (2015). Revisiting approximate reanalysis in topology optimization: 
on the advantages of recycled preconditioning in a minimum weight procedure. 
Structural and Multidisciplinary Optimization, 51(1), 41-57.

The codes are now freely available, please reference the relevant publications in your research work.

Three types of codes:
  1) top3d*.m = 3-D minimum compliance topopt using MGCG solver:
  
		1.1) basic
		1.2) with sensitivity monitoring
		
		
		
  
  2) minV*.m = 2-D minimum volume topopt:
     a) basic
     b) reanalysis
     c) with matrix-free reanalysis
  3) minW3d*.m = 3-D minimum volume topopt using MGCG solver:
     a) basic
     b) matrix-free with recycling preconditioner
     c) matrix-free with recycling preconditioner and sensitivity monitoring

Extensions and improvements are strongly encouraged!!!
