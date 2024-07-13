# SATINN
Software for Analysis of Testis Images with Neural Networks

This repository hosts the code and UI as described for SATINN 1 in [Yang et al., 2022](https://academic.oup.com/bioinformatics/article/38/23/5288/6754803) and for SATINN 2 in Yang et al., 2024 (manuscript in preparation).

&nbsp;

## Download SATINN

Navigate to the `deploy` folder in this repository to download the Matlab-based app.

* [Recommended] Users **with Matlab R2023b or newer** already installed on their device should download the contents of `for_redistribution_files_only`. No additional install is required.
* Users **without Matlab** should download the installer from `for_redistribution` and will need to install Matlab components.
* Users **with an earlier Matlab release** should consider updating to at least R2023b to ensure future compatibility, then follow the Recommended step above.

Additionally, **all users** should download our most recent [demo data from Figshare](https://figshare.com/articles/dataset/SATINN_-_demo_files_and_neural_networks_r_230516_/22853645), which contains a sample test image as well as our pre-trained neural networks. In this RAR, the `datasets` and `neuralnets` folders should be extracted in the same directory as `SATINN.exe` as shown below.

![Example_Install](https://github.com/conradlab/SATINN/assets/43147040/8061e31a-2e24-434e-826f-afa4a9fd1590)

## Demo instructions

Once SATINN is properly installed, running `SATINN.exe` should bring up the Inputs window that looks like this:

![SATINN_Inputs](https://github.com/conradlab/SATINN/assets/43147040/83b6a70d-f68a-4d23-958f-0af832ed5695)

To try out our data:

* Raw image - select `MS124RS101_med.tif` (from the `datasets` folder)
* Cell segmentation - select `MS124RS101_cpm_filt_med.tif`
* Image type - This a brightfield image fixed with PAS, so choose 'Brightfield (PAS)'.
* For our data, you can leave the Resolution options as is.
* Cell type neural network - select `neuralnets/bf_cnetr-230408.mat`. (Ensure `neuralnets` is in the same directory as `SATINN.exe` as shown at the end of the [Download](#download-satinn) section.)
* Tubule stage neural network - select `neuralnets/bf_tnetr-230314.mat`

Then click Process. The rest is automated, and will take a few minutes. (On a 32 GB RAM computer, our dataset takes about 3 minutes to fully process.) Once the Process light turns green, you can view the results in the [Cell Outputs] and [Tubule Outputs] panels, as well as save the full datatables and more using the [Save...] option.

![SATINN_Outputs](https://github.com/conradlab/SATINN/assets/43147040/d7804432-d470-4322-9285-9a121b411a56)

Note: The SATINN App is still in development. Please report any app-breaking bugs by opening a New Issue.

&nbsp;

## A note on cell segmentation

If you wish to use SATINN on your own data, we strongly recommend augmenting your input by using Cellpose to perform cell segmentation. For the Yang et al. 2022 manuscript, we used the command line version of Cellpose to process all of our images. This command uses the settings (usually default) that allow exact reproduction:

     python -m cellpose --dir <<your dir>> --pretrained_model cyto --chan 0
     --save_tif --no_npy --use_gpu --batch_size 8 --flow_threshold 0

For more info on Cellpose, including instructions for installation, see: http://www.cellpose.org.
