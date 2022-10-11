# SATINN
Software for Analysis of Testis Images with Neural Networks

This repository hosts code and UI for SATINN in Yang et al., 2022.

See the `deploy` folder to download the Matlab-based app. A Matlab subscription is not required, though Matlab itself might still need to be downloaded.

The SATINN UI contains a built-in cell classification neural network. However, the tubule classification netowrk is hosted separately and can be accessed using the Figshare link below. To use it in the SATINN app, in the Tubule stages pre-trained neural network drop-down, choose Custom and specify this file.

Tubule neural network Matlab file: https://figshare.com/articles/dataset/Mouse_Tubule_Classification_Neural_Network_File/21313107

A demo image dataset can be download from Figshare below, which can be used directly with either the app or the demo code. A subset of the dataset is also provided to reduce load times when testing the functionality of SATINN:

Full: https://figshare.com/articles/dataset/MS36R1_SEC2B_rar/19619010

Subsection: https://figshare.com/articles/dataset/MS36R1_SEC2B_demo_rar/19822441



&nbsp;

The MATLAB script demo.m is available as an alternative to the UI and contains all steps necessary to reproduce the results in Yang et al. 2022. Although the results in the manuscript are a culmination of multiple image sources, demo.m uses the image set available above, MS36R1_SEC2B, which contributed to neural network training and validation, and wildtype cell and tubule analysis in the paper.

If using the raw demo script, check that you have the following folders and files:
 * A `datasets` folder containing the `MS36R1_SEC2B` dataset (see above). This folder should then contain four images: Hoe, Hoe_cp_masks, Acta2, and Acrv1
 * A `functions` folder containing 7 supplementary functions for use with this demo.
 * The file `cellnet_h-XXXXXX.mat` or `cnetr-XXXXXX.mat`, which contains the cell classification neural net.
 
 Be sure to add the `functions` folder to path, or move the functions there to your current path. You will also need to place the downloaded image data in a folder called `datasets`. For processing additional data, use the same format as the files shown.

 For the manuscript, we used the command line version of Cellpose to segment all of our cell images. This command uses the settings (usually default) that allow exact reproduction:

     python -m cellpose --dir <<your dir>> --pretrained_model cyto --chan 0
     --save_tif --no_npy --use_gpu --batch_size 4

 However, as Cellpose does take a while to execute on large images (ours is on the order of 25000x25000), we have included the output of Cellpose as part of the demo dataset, for convenience.
  
 For more info on Cellpose, including instructions for installation, see: http://www.cellpose.org.
