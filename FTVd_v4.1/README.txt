************************************************************************
            FTVd_v4: Fast Total Variation Deconvolution (version 4)
************************************************************************

Copyright (C) 2009 Junfeng Yang, Yin Zhang, Wotao Yin and Yilun Wang

1). Get Started
===================

   Run any demo code listed below:

         demoTVL2  : Demo TV/L2 solve
         demoMTVL2 : Demo multichannel TV/L2 solve
         demoTVL1  : Demo TV/L1 solve
         demoMTVL1 : Demo multichannel TV/L1 solve

   or run demoAll to demostrate all.

2). Introduction
====================
   
FTVd_v4 keeps all the features of FTVd_v3.0 and is a major upgrade.
In this version, we implemented the Alternating Direction Method (ADM),
which dates back to Glowinski and Marrocco (1975) and Gabay and Mercier
(1976). 

Compared with ver 3.0, this version avoids the ill-conditioning caussed by 
large penalty paramter by considering the augmented Lagrangian and ADM. 
As a result, FTVd_v4 has a much simply parameter selection rule and 
behaves very efficiently and roubustly.

For more introductory material about FTVd, see READ_FTVd_v3.0, which is a 
copy of README in FTVd_v3.0.
 

3). Usage
====================

FTVd_v4 is called in the following way:

               out = FTVd_v4(Bn,H,mu,str), 

where 

     Bn  -- a blurry and noiy observation,
     H   -- a convolution kernel, 
     mu  -- Lagrangian multiplier (positive number)
     str -- a string used to specify a model. 
            * str = 'L2' or 'l2' to specify a L2 model (for white noise),
            * str = 'L1' or 'l1' to specify a L1 model (for impulsive noise).  

Users can also specify algorithm paramters and call FTVd_v4 like

               out = FTVd_v4(Bn,H,mu,str,opts), 

where opts is a structure with some fields. For example, when L1 model
is to be solved, users may specify fields listed below:
  
% opts -- a structure containing algorithm parameters {default}
        * opst.beta1   : a positive constant {5}
        * opst.beta2   : a positive constant {20}
        * opst.gamma   : a constant in (0,1.618] {1.618}
        * opst.maxitr  : maximum iteration number {500}
        * opst.relchg  : a small positive parameter which controls
                         stopping rule of the code. When the
                         relative change of X is less than
                         opts.relchg, then the code stops. {1.e-3}
        * opts.print   : print inter results or not {0}
For L2 models, .beta1 and .beta2 should be replaced by one field .beta. For 
more explanation about FTVd, see each core solver in folder "coresolver"
and also check READ_FTVd_v3.0.



4). References
====================

    For algorithmic details, such as continution on penalty parameters and 
 optimality conditions, see references:

  [1] Y. Wang, J. Yang, W. Yin and Y. Zhang, "A New Alternating Minimiza-
      tion Algorithm for Total Variation Image Reconstruction", 
      SIAM Journal on Imaging Sciences, 1(3), 248-272, 2008.

  [2] J. Yang, W. Yin, Y. Zhang and Y. Wang, "A Fast Algorithm for Edge-
      Preserving Variational Multichannel Image Restoration", SIAM Journal 
      on Imaging Sciences, 2(2), 569-592, 2009.

  [3] J. Yang, Y. Zhang and W. Yin, "An efficient TVL1 algorithm for 
      deblurring multichannel images corrupted by impulsive noise",
      SIAM Journal on Scientific Computing, 31(4), 2842-2865, 2009.

  [4] M. Tao, J. Yang and B. He, "Alternating direction algorithms for total 
      variation deconvolution in image reconstruction ", TR09-18, Department 
      of Mathmatics, Nanjing University, August, 2009, Available at Optimization
      online: http://www.optimization-online.org/DB_HTML/2009/11/2463.html
     

5). Contact Information
=======================

FTVd is available at: http://www.caam.rice.edu/~optimization/L1/ftvd/

Please feel free to e-mail the following authors with any comments 
or suggestions:

Junfeng Yang,  Depart. Math., Nanjing Univ.,  <jfyang@nju.edu.cn>


6).  Copyright Notice
====================

   FTVd is free software; you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free 
Software Foundation; either version 3 of the License, or (at your option) 
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details at
<http://www.gnu.org/licenses/>. 

