function setup_tubclass2(im, ann, masks, outdir)

% % example:
% im = cat(3, v, h, a);
% setup_tubclass2(im, ...
%     'annotations.csv', ...
%     'tub_ROI.tif', ...
%     '../outputs/outdir');

fprintf('Read annotations\n')
imann = readtable(ann);

fprintf('Import masks\n')
imabwld = imread(masks);

fprintf('Calculate regionprops\n')
imaR = regionprops(imabwld);

fprintf('Extract ROIs\n')
extract_objects(im, imabwld, imaR, imann, 2000, 'outdir', outdir);

fprintf('Adjust image size\n')
f = dir([outdir, '*.tif']);
for i = 1:size(f, 1)
    ff = fullfile(f(i).folder, f(i).name);
    ffi = imread(ff);
    imwrite(imresize(ffi, 0.25), ff);
end