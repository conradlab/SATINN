function setup_cellclass(imp, imh_cp_masks, imabwld, indices, varargin)

% function setup_cellclass(hfile, afile, abwldfile, indices, varargin)
% 
% % extract cells for annotation from specific tubules
% % example:
% % hfile = 'C:\Users\yangra\Documents\R\projects\ima\ana\datasets\MS36R1_SEC2B\MS36R1_SEC2B-HOE.tif';
% % afile = 'C:\Users\yangra\Documents\R\projects\ima\ana\datasets\MS36R1_SEC2B\MS36R1_SEC2B-ACTA2.tif';
% % abwldfile = 'C:\Users\yangra\Documents\MATLAB\outputs\groundtruth\refs\MS36R1_SEC2B_tub_masks.tif';
% % annos = readtable('C:\Users\yangra\Documents\R\projects\ima\ana\datasets\MS36R1_SEC2B\MS36R1_SEC2B_tub_stages.csv');
% % tubuleinds = randsample(annos.Tubule(cellfun(@(x) ismember(x, {'1'}), annos.Stage)), 3); % these are 3 random samples of stage 1
% % setup_cellclass(hfile, afile, abwldfile, tubuleinds, 'outdir', '../outputs');
% 
% p = inputParser;
% addRequired(p, 'hfile', @(x) isstring(x) | ischar(x));
% addRequired(p, 'afile', @(x) isstring(x) | ischar(x));
% addRequired(p, 'abwldfile', @(x) isstring(x) | ischar(x));
% addRequired(p, 'indices', @(x) isnumeric(x));
% addParameter(p, 'outdir', [], @(x) isstring(x) | ischar(x));

p = inputParser;
addRequired(p, 'imp');
addRequired(p, 'imh_cp_masks');
addRequired(p, 'imabwld');
addRequired(p, 'indices');
addRequired(p, 'outdir', @(x) isstring(x) | ischar(x));
addRequired(p, 'prefix', @(x) isstring(x) | ischar(x));

parse(p, imp, imh_cp_masks, imabwld, indices, varargin{:});
im = p.Results.imp;
imhbwl = p.Results.imh_cp_masks;
imabwld = p.Results.imabwld;
indices = p.Results.indices;
outdir = p.Results.outdir;
name = p.Results.prefix;

% fprintf(' Debug: parsed outdir as: %s\n', outdir);
% fprintf(' Debug: parsed prefix as: %s\n', name);

% 
% parse(p, hfile, afile, abwldfile, indices, varargin{:});
% hfile = p.Results.hfile;
% afile = p.Results.afile;
% abwldfile = p.Results.abwldfile;
% indices = p.Results.indices;
% outdir = p.Results.outdir;
% 
% if isempty(outdir)
%     outdir = '../outputs';
% end
% 
% imh = imread(hfile);
% ima = imread(afile);
% 
% [fp, ~, ~] = fileparts(hfile);
% % fp2 = strsplit(fp, '\');
% fp2 = regexp(fp, '[\/|\\]', 'split');
% name = fp2{end};
% 
% imhp = imadjust(imtophat(imh, strel('disk', 100)));
% imap = imadjust(imtophat(ima, strel('disk', 100)));
% 
% imabwld = imread(abwldfile);

% fprintf(' Image processing\n')
imaR = regionprops(imabwld);

% construct internal bwdist and rescale
imbwdist = bwdist(~imabwld);
for i = 1:size(imaR, 1)
    imbwdist(imabwld == i) = rescale(imbwdist(imabwld == i), 0, 65535);
end

% construct rgb image of relevant info
% im2 = cat(3, imp(:, :, 2), imbwdist, imp(:, :, 3));
im = cat(3, im, imbwdist);

% isolate ROIs of interest
tmp = reshape([imaR.BoundingBox], [4, size(imaR, 1)])';
pval = round(max(tmp(:, 3:4), [], 'all'));

im2Out = cell(length(indices), 1);
imaOut = cell(length(indices), 1);
imhOut = cell(length(indices), 1);
im = padarray(im, [pval, pval], 0, 'both');
imabwldp = padarray(imabwld, [pval pval], 0, 'both');
imhbwlp = padarray(imhbwl, [pval pval], 0, 'both');
fprintf(' Extract single tubules    0%%');
if ~exist(fullfile(outdir, name), 'dir')
    mkdir(fullfile(outdir, name))
end
for i = 1:length(indices)
    fprintf('\b\b\b\b\b% 4.0f%%', i/length(indices)*100);
    x = indices(i);
    if ~exist(fullfile(outdir, [name, '-', sprintf('%03.0f', x)]), 'dir')
        mkdir(fullfile(outdir, [name, '-', sprintf('%03.0f', x)]));
    end
    
    im2Out{i, 1} = im(ceil(imaR(x).BoundingBox(2)+pval-50):ceil(imaR(x).BoundingBox(2)+pval+imaR(x).BoundingBox(4)+50), ...
                       ceil(imaR(x).BoundingBox(1)+pval-50):ceil(imaR(x).BoundingBox(1)+pval+imaR(x).BoundingBox(3)+50), :);
    imwrite(im2Out{i, 1}(:, :, [2 4 3]), fullfile(outdir, [name, '-', sprintf('%03.0f', x)], [name, '-', sprintf('%03.0f', x), '_RGB.tif']));
    imwrite(im2Out{i, 1}(:, :, 2), fullfile(outdir, [name, '-', sprintf('%03.0f', x)], [name, '-', sprintf('%03.0f', x), '_HOE.tif']));
    
    imaOut{i, 1} = imabwldp(ceil(imaR(x).BoundingBox(2)+pval-50):ceil(imaR(x).BoundingBox(2)+pval+imaR(x).BoundingBox(4)+50), ...
                            ceil(imaR(x).BoundingBox(1)+pval-50):ceil(imaR(x).BoundingBox(1)+pval+imaR(x).BoundingBox(3)+50), :);
    imaOut{i, 1}(imaOut{i, 1} ~= x) = 0;
    imaOut{i, 1} = logical(imaOut{i, 1});
    imwrite(imaOut{i, 1}, fullfile(outdir, [name, '-', sprintf('%03.0f', x)], [name, '-', sprintf('%03.0f', x), '_tub_mask.tif']));
    
    tmp = im2Out{i, 1}(:, :, 1:3);
    tmp(repmat(~imaOut{i, 1}, 1, 1, 3)) = 0;
    imwrite(tmp, fullfile(outdir, name, [name, '-', sprintf('%03.0f', x), '.tif']));
    
    % isolate cells within tubules
    imhOut{i, 1} = imhbwlp(ceil(imaR(x).BoundingBox(2)+pval-50):ceil(imaR(x).BoundingBox(2)+pval+imaR(x).BoundingBox(4)+50), ...
                           ceil(imaR(x).BoundingBox(1)+pval-50):ceil(imaR(x).BoundingBox(1)+pval+imaR(x).BoundingBox(3)+50), :);
    cellconv = nonzeros(unique(imhOut{i, 1}));
    for j = 1:size(cellconv, 1)
        imhOut{i, 1}(imhOut{i, 1} == cellconv(j)) = j;
    end
    imwrite(uint16(imhOut{i, 1}), fullfile(outdir, [name, '-', sprintf('%03.0f', x)], [name, '-', sprintf('%03.0f', x), '_HOE_cp_masks.tif']));
    cellconv = array2table([(1:size(cellconv, 1))', cellconv]);
    cellconv.Properties.VariableNames = {'LocalID', 'GlobalID'};
    writetable(cellconv, fullfile(outdir, [name, '-', sprintf('%03.0f', x)], [name, '-', sprintf('%03.0f', x), '_HOE_cellID_conversions.csv']));
end
fprintf('\n');