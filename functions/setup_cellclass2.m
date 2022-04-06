function [imR2, centroids] = setup_cellclass2(imfile, imbwlfile, imabwlfile, varargin)

% setup more cell class stuff
% example
% imfile = '../savefiles/cellclass/MS36R1_SEC2B-006/MS36R1_SEC2B-006_HOE.tif';
% imbwlfile = '../savefiles/cellclass/MS36R1_SEC2B-006/MS36R1_SEC2B-006_HOE_cp_masks.tif';
% imabwlfile = '../savefiles/cellclass/MS36R1_SEC2B-006/MS36R1_SEC2B-006_tub_mask.tif';
% centroids = setup_cellclass2(imfile, imbwlfile, imabwlfile);
% writematrix(centroids, '../savefiles/cellclass/MS36R1_SEC2B-006/MS36R1_SEC2B-006_centroids.csv');

p = inputParser;
addRequired(p, 'imfile', @(x) isstring(x) | ischar(x));
addRequired(p, 'imbwlfile', @(x) isstring(x) | ischar(x));
addRequired(p, 'imabwlfile', @(x) isstring(x) | ischar(x));
addParameter(p, 'savefiles', true, @(x) islogical(x));
addParameter(p, 'cellconv', [], @(x) isempty(x) | isstring(x) | ischar(x));

parse(p, imfile, imbwlfile, imabwlfile, varargin{:});
imfile = p.Results.imfile;
imbwlfile = p.Results.imbwlfile;
imabwlfile = p.Results.imabwlfile;
savefiles = p.Results.savefiles;
cellconv = p.Results.cellconv;

if ~isempty(cellconv)
    cellconv = readtable(cellconv);
end

[fp, ~, ~] = fileparts(imfile);
fp2 = regexp(fp, '[\/|\\]', 'split');
name = fp2{end};

im = imread(imfile);
imbwl = imread(imbwlfile);
imabwl = imread(imabwlfile);

imR = regionprops(imbwl);

tmp = reshape([imR.Centroid], [2, size(imR, 1)])';

% keep cells only within tubule bounds (logical indexing)
filter = imabwl(sub2ind(size(imabwl), round(tmp(:, 2)), round(tmp(:, 1))));
imR2 = imR(filter);
if ~isempty(cellconv)
    cellconv.isintubule = filter;
end
centroids = reshape([imR2.Centroid], [2, size(imR2, 1)])';
if savefiles
    writematrix(centroids, fullfile(fp, [name, '_centroids.csv']));
    if ~isempty(cellconv)
        writetable(cellconv, fullfile(fp, [name, '_HOE_cellID_conversions.csv']));
    end
end

f = figure;
imshow(im)
hold on
scatter(centroids(:, 1), centroids(:, 2), 'filled');
text(centroids(:, 1), centroids(:, 2), ...
    cellfun(@(x) num2str(x), num2cell(find(imabwl(sub2ind(size(imabwl), round(tmp(:, 2)), round(tmp(:, 1)))))), 'uniformoutput', false), ...
    'color', 'r', 'fontname', 'fixedwidth', 'fontsize', 10, 'fontweight', 'bold', 'horizontalalignment', 'center');
title(name, 'interpreter', 'none');
if savefiles
    saveas(f, fullfile(fp, [name, '_centroids_overlay.tif'])); 
    close(f)
end
end