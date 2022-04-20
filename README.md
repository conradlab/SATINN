# SATINN
Software for Analysis of Testis Images with Neural Networks

This repository hosts code for SATINN. The MATLAB code demo.m contains all 
steps necessary to reproduce the results in Yang et al. 2022. Although the results 
in the manuscript are a culmination of multiple image sources, demo.m
uses a single image set, MS36R1_SEC2B, which contributed to neural
network training and validation, and wildtype cell and tubule analysis in
the paper.

The dataset can be download from Figshare at the following address:
https://figshare.com/articles/dataset/MS36R1_SEC2B_rar/19619010

Before you begin, check that you have the following folders and files:
 * A 'datasets' folder containing the 'MS36R1_SEC2B' dataset (see above)
  This folder should then contain four images: Hoe, Hoe_cp_masks, Acta2,
   and Acrv1
 * A 'functions' folder containing 7 supplementary functions for use with
   this demo.
 * The file 'cellnet_h-XXXXXX.mat', which contains the latest neural net
   created for cell classification.
 
 Be sure to add the `functions` folder to path, or move the functions
 there to your current path.

 You also need to place your dataset folder (e.g. `MS36R1_SEC2B`) in the
 folder called `datasets`. For processing additional data, use the same
 format as the files shown.

 For the manuscript, we used the command line version of cellpose to
 segment all of our cell images, including the one used in demo.m. This
 includes the settings (usually default) that allow exact reproduction:

     python -m cellpose --dir <<your dir>> --pretrained_model cyto --chan 0
     --save_tif --no_npy --use_gpu --batch_size 4

 However, as Cellpose does take a while to execute on large images (ours
 is approximately 25000x25000), we have included the output of Cellpose
 as part of the demo dataset, for convenience.
  
 For more info on Cellpose, including instructions for installation, see:
 http://www.cellpose.org.


