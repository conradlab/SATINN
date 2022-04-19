# SATINN
Software for Analysis of Testis Images with Neural Networks

This repository hosts code for SATINN. The MATLAB code demo.m contains all 
steps necessary to reproduce the results in Yang et al. 2022. Although the results 
in the manuscript are a culmination of multiple image sources, demo.m
uses a single image set, MS36R1_SEC2B, which contributed to neural
network training and validation, and wildtype cell and tubule analysis in
the paper.

Before you begin, check that you have the following folders and files:
 * A 'datasets' folder containing the 'MS36R1_SEC2B' dataset.
  This folder should then contain four images: Hoe, Hoe_cp_masks, Acta2,
   and Acrv1
 * A 'functions' folder containing 7 supplementary functions for use with
   this demo.
 * The file 'cellnet_h-XXXXXX.mat', which contains the latest neural net
   created for cell classification.
 
 Be sure to add the `functions` folder to path, or move the functions
 there to your current path.

 Datasets in the `dataset` folder should look like this, in case you want
 to analyze additional datasets:

 Working directory (specify this in the variable `wd`, e.g. '../demo/datasets')
  + Dataset1
      Dataset1-HOE.tif
      Dataset1-ACTA2.tif
      Dataset1-ACRV1.tif
      Dataset1-HOE_cp_masks.tif (Hoechst masks from Cellpose or another
                                 segmentation output of your choice)
  + Dataset2
      etc.

 For the manuscript, we used the command line version of cellpose to
 segment all of our cell images, including the one used in demo.m. This
 includes the settings (usually default) that allow exact reproduction:

   python -m cellpose --dir <<your dir>> --pretrained_model cyto --chan 0
   --save_tif --no_npy --use_gpu --batch_size 4

 For more info on Cellpose, including instructions for installation, see:
 http://www.cellpose.org.


