% TICAT: Testis Image Classification and Analysis Tool
%
% This is a demo code for TICAT. It contains all steps necessary to
% reproduce the results in Yang et al. 2022. Although the results discussed
% in the manuscript are a culmination of multiple image sources, this demo
% uses a single image set, MS36R1_SEC2B, which contributed to neural
% network training and validation, and wildtype cell and tubule analysis in
% the paper.
%
% Before you begin, check that you have the following folders and files:
% * A 'datasets' folder containing the 'MS36R1_SEC2B' dataset.
%    This folder should then contain four images: Hoe, Hoe_cp_masks, Acta2,
%    and Acrv1
% * A 'functions' folder containing 7 supplementary functions for use with
%   this demo.
% * The file 'cellnet_h-XXXXXX.mat', which contains the latest neural net
%   created for cell classification.
% 
% Be sure to add the `functions` folder to path, or move the functions
% there to your current path.
%
% Datasets in the `dataset` folder should look like this, in case you want
% to analyze additional datasets:
%
% Working directory (specify this in the variable `wd`, e.g. '../demo/datasets')
%  + Dataset1
%      Dataset1-HOE.tif
%      Dataset1-ACTA2.tif
%      Dataset1-ACRV1.tif
%      Dataset1-HOE_cp_masks.tif (Hoechst masks from Cellpose or another
%                                 segmentation output of your choice)
%  + Dataset2
%      etc.
%
% For this manuscript, we used the command line version of cellpose to
% segment all of our cell images, including the one used in this demo. This
% includes the settings (usually default) that allow exact reproduction:
%
%   python -m cellpose --dir <<your dir>> --pretrained_model cyto --chan 0
%   --save_tif --no_npy --use_gpu --batch_size 4
%
% For more info on Cellpose, including instructions for installation, see:
% http://www.cellpose.org.

%% Neural network training -- cell classification
%
% You can also use the pre-trained network included in this demo, located
% in the file: ../demo/cellnet_h-xxxxxx.mat (where xxx is the date of the
% most recent version). If you do, you can skip this section.
%
% This section specifically covers cell type classification, but you can
% train a NN to classify any kind of image data this way. The architecture
% might need to be adjusted based on the expected features, but the overall
% procedure remains the same.

% Load image data (you can modify this or input your own annotated data
% here).

% This is the fraction of the annotated data that will be used for
% training. The remainder will be used to validate the result.
pctTrainFiles = 0.7;
[dsTrain, dsValidation] = splitEachLabel(ds, pctTrainFiles, 'randomize');

% Define network architecture. For more information on designing and
% training neural networks using Matlab's Deep Learning Toolbox, see their
% resource: https://www.mathworks.com/products/deep-learning.html
layers = [
    imageInputLayer(size(imread(ds.Files{1, 1})))
    
    convolution2dLayer(2,8,'Padding','same', 'stride', 2)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(2,16,'Padding','same', 'stride', 2)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(2,32,'Padding','same', 'stride', 2)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(length(unique(ds.Labels)))
    softmaxLayer
    classificationLayer];

% Set training conditions
options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',20, ...
    'Shuffle','every-epoch', ...
    'ValidationData',dsValidation, ...
    'ValidationFrequency',3, ...
    'ValidationPatience',40, ...
    'Verbose',false, ...
    'Plots','training-progress');

% Perform NN training
net = trainNetwork(dsTrain, layers, options);

% Evaluate trained network
cnnValidationFigs(net, dsValidation);

%% Image processing

% Load directory containing datasets. Make sure the minimum files are
% present (see instruction header). All additional files generated here
% will overwrite existing files with the same name.
wd = '../demo/datasets';
f = dir(wd);
f = f([f.isdir]);

% If you don't want to process all files here, specify them by their
% indices in `f`. The default range for all files is 3:size(f, 1).
for F = 3:size(f, 1)
    % Print dataset status update
    c = clock; fprintf('%02.0f:%02.0f:%02.0f: Working on %s...\n', c(4), c(5), round(c(6)), f(F).name);
    
    % Read raw image files
    imh = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-HOE.tif']));
    ima = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-ACTA2.tif']));
    imv = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-ACRV1.tif']));
    
    % Initial tubule processing. Normalization followed by tubule
    % segmentation and image assembly.
    fprintf('  Process tubules\n');
    [im, imabwl, ~] = setup_tubclass(imh, ima, imv);
    
    % Occasionally free up memory to avoid out of memory errors for large 
    % datasets when running locally
    clear imh ima imv
    
    % Assemble tubule statistics. Dilation strel accounts for the Acta2
    % dilation factor necessary to get good segmentation (basically just
    % undoes a processing step from `setup_tubclass`).
    fprintf('  Label tubules\n');
    imabwld = imdilate(bwlabel(imabwl), strel('disk', 30));
    clear imabwl
    % Save tubule segmentation mask info
    imwrite(uint8(imabwld), fullfile(f(F).folder, f(F).name, [f(F).name, '_tub_masks.tif']));
    
    % Initial cell processing block
    fprintf('  Process + extract all cells for classification\n');
    % Generate relative apical-basal distance data for each tubule. This
    % gives us a normalized method to locate cells within and between
    % tubules.
    imbwdist = bwdist(~imabwld);
    imaR = regionprops(imabwld);
    for i = 1:size(imaR, 1)
        imbwdist(imabwld == i) = rescale(imbwdist(imabwld == i), 0, 65535);
    end
    % Establish cell image containing Hoe, A/B position, and Acta2 data
    % (preserved from previous builds but currently unused).
    im2 = cat(3, im(:, :, 2), imbwdist, im(:, :, 3));
    % Read cellpose segmentation data and build stats table
    imhbwl = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-HOE_cp_masks.tif']));
    imhR = regionprops(imhbwl);
    % Extract individual cell images from entire slide (warning: this might
    % take a while if there are a lot of cells).
    ims = extract_objects(im2, [], imhR, [], 50);
    % Generate and save cell image stack as Matlab file. Specifically, the
    % `cells` variable will be a 50x50x3xN stack of N cell images. To
    % access this data later, call `load('<dataset_name>_cells.mat')` which
    % will read in `cells` to your workspace.
    cells = cat(4, ims{:, 1});
    save(fullfile(f(F).folder, f(F).name, [f(F).name, '_cells.mat']), 'cells', '-v7.3')
    
    % Cleanup again to avoid out of memory errors for large datasets.
    clear im2 imbwdist ims cells
    
    
    % Create folders for each individual tubule found and analyzed within
    % this particular slide. Some of this info will be used in the next
    % section to build cell stats.
    fprintf('  Create individual tubule folders\n');
    setup_cellclass(im, imhbwl, imabwld, (1:size(imaR, 1))', 'prefix', f(F).name, 'outdir', fullfile(f(F).folder, f(F).name, 'tubules'));
    % zip tubule images for annotations
    if isfolder(fullfile(f(F).folder, f(F).name, 'tubules', f(F).name))
        zip(fullfile(f(F).folder, f(F).name, 'tubules', [f(F).name, '.zip']), fullfile(f(F).folder, f(F).name, 'tubules', f(F).name));
    end
    delete(fullfile(f(F).folder, f(F).name, 'tubules', f(F).name, '*'));
    rmdir(fullfile(f(F).folder, f(F).name, 'tubules', f(F).name));
    h = dir(fullfile(f(F).folder, f(F).name, 'tubules'));
    h = h([h.isdir]);
    for H = 3:size(h, 1)  % set this index to 4 if not zipping/removing tubule annotation folder
        loc = fullfile(h(H).folder, h(H).name);
        setup_cellclass2(fullfile(loc, [h(H).name, '_HOE.tif']), ...
            fullfile(loc, [h(H).name, '_HOE_cp_masks.tif']), ...
            fullfile(loc, [h(H).name, '_tub_mask.tif']), ...
            'cellconv', fullfile(loc, [h(H).name, '_HOE_cellID_conversions.csv']));
    end
    clearvars -except f wd
end

%% Classification

% Load pre-trained neural net for cells, as `net`. You can also skip this
% line if you have trained your own neural net.
load('../demo/cellnet_h-210731.mat');

% Same as above, specifies directory containing data to be analysed
wd = '../demo/datasets';
f = dir(wd);
f = f([f.isdir]);

% If you don't want to process all files here, specify them by their
% indices in `f`. The default range for all files is 3:size(f, 1).
for F = 3:size(f, 1)
    % Print dataset status update
    c = clock; fprintf('%02.0f:%02.0f:%02.0f: Working on %s\n', c(4), c(5), round(c(6)), f(F).name);
    
    % Load stack of individual cells and classify using pre-trained neural
    % net
    load(fullfile(f(F).folder, f(F).name, [f(F).name, '_cells.mat']));
    c = clock; fprintf('  %02.0f:%02.0f:%02.0f: Classify\n', c(4), c(5), round(c(6)));
    [p, s] = classify(net, cells(:, :, 1, :));

    % Build cell stats table
    c = clock; fprintf('  %02.0f:%02.0f:%02.0f: Building cell data table\n', c(4), c(5), round(c(6)));
    TC = array2table(s);
    TC.Properties.VariableNames = cellstr(net.Layers(15, 1).Classes);
    TC.GlobalID = (1:size(s, 1))'; % unique indicator of cells w.r.t. whole image 
    TC.LocalID = zeros(size(s, 1), 1); % indicator of cells w.r.t. contained tubule (0 if not in tubule)
    TC.InTubule = zeros(size(s, 1), 1); % indicator of tubule that contains this particular cell (0 if not in tubule)
    g = dir(fullfile(f(F).folder, f(F).name, 'tubules'));
    g = g([g.isdir]);
    for G = 4:size(g, 1)
        cnv = readtable(fullfile(g(G).folder, g(G).name, [g(G).name, '_HOE_cellID_conversions.csv']));
        t = str2double(g(G).name(end-2:end));
        cnv = cnv(logical(cnv.isintubule), :);
        for i = 1:size(cnv, 1)
            TC.LocalID(TC.GlobalID == cnv.GlobalID(i)) = cnv.LocalID(i);
            TC.InTubule(TC.GlobalID == cnv.GlobalID(i)) = t;
        end
    end
    TC = movevars(TC, {'GlobalID', 'LocalID', 'InTubule'}, 'before', 1);
    TC.Celltype = p; % predicted cell type
    TC.Celltype = categorical(TC.Celltype);
    TC.RelPos = zeros(size(s, 1), 1); % A/B position
    for i = 1:size(s, 1)
        tmp = squeeze(cells(:, :, :, i));
        TC.RelPos(i) = mean(tmp(:, :, 2), 'all');
    end
    
    % Incorporate physical cell stats into table
    imhbwl = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-HOE_cp_masks.tif']));
    imhR = regionprops(imhbwl, {'Area', 'BoundingBox', 'Centroid', 'Circularity', 'Eccentricity', ...
        'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});
    imhT = struct2table(imhR);
    TC = [TC, imhT]; %#ok<AGROW>
    
    % Calculate cell relative orientation
    imabwld = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '_tub_masks.tif']));
    imabwldR = regionprops(imabwld);
    TC.RelOrient = zeros(size(s, 1), 1); % cell orientation relative to A/B
    for i = 1:size(s, 1)
        if TC.InTubule(i) ~= 0
            v = TC.Centroid(i, :) - imabwldR(TC.InTubule(i)).Centroid;
            TC.RelOrient(i) = min([abs(TC.Orientation(i) - atan2(v(2), v(1))), abs(TC.Orientation(i) - atan2(-v(2), -v(1)))]);
        end
    end
    
    
    % Save cell data table
    writetable(TC, fullfile(f(F).folder, f(F).name, [f(F).name, '_cellstats_h.csv']));
    clearvars -except f net wd
end

%% Evaluate classification

% We can also directly visualize some of the cell classifier's calls. Here
% is an example of one way to do so, but there are many ways to make
% various figures from the stats table (..._cellstats_h.csv).

% Same as above, specifies directory containing data to be analysed
wd = '../demo/datasets';
f = dir(wd);
f = f([f.isdir]);

% If you don't want to process all files here, specify them by their
% indices in `f`. The default range for all files is 3:size(f, 1).
for F = 3:size(f, 1)
    % Import necessary information
    TC = readtable(fullfile(f(F).folder, f(F).name, [f(F).name, '_cellstats_h.csv']));
    imh = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-HOE.tif']));
    imhbwl = imread(fullfile(f(F).folder, f(F).name, [f(F).name, '-HOE_cp_masks.tif']));
    TC.Celltype = categorical(TC.Celltype);

    % Show cell segmentation output over original image
    L = label2rgb(imhbwl, 'parula', 'k', 'shuffle');
    figure
    imshow(imh)
    hold on
    h = imshow(L);
    h.AlphaData = 0.4;
    title(sprintf('%s Segmentation', f(F).name));
    
    % Show classified cell type over original image
    Lnum = uint8(zeros(size(imhbwl)));
    Lcats = categories(TC.Celltype);
    for i = 1:size(Lcats, 1)
        val = TC.GlobalID(TC.Celltype == Lcats{i});
        Lnum(ismember(imhbwl, val)) = i;
    end
    L = label2rgb(Lnum, 'parula', 'k');
    figure
    imshow(imh)
    hold on
    h = imshow(L);
    h.AlphaData = 0.6;
    title(sprintf('%s Classification', f(F).name));
end