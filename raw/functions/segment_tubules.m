function [IMbwl, IMB, IMR] = segment_tubules(IM, varargin)

%% Segment tubule (Acta2)


% Inputs:
% IM: Image to be segmented

% Optional inputs, specified as name-value pair arguments:
% d_val: dilation value, initial image dilation value in attempt to join
% boundaries that are not fully closed. Default value 25 for 20k images.
% o_val: open value, to remove small objects and holes in boundaries,
% default value 100 for 20k images.
% figures: Option to display figures for validation purposes. Default
% false.
% verbose: Adds additional information during function runtime. Default
% false.

% Example:
% IM = imread(...);
% [IMbwl, IMB, IMR] = segment_tubules(IM, 'd_val', 25, 'o_val', 100, 'figures', true, 'verbose', true);


p = inputParser;
addRequired(p, 'IM');
addParameter(p, 'd_val', 25, @(x) isnumeric(x));
addParameter(p, 'o_val', 100, @(x) isnumeric(x));
addParameter(p, 'figures', false, @(x) islogical(x));
addParameter(p, 'verbose', false, @(x) islogical(x));

parse(p, IM, varargin{:});
IM = p.Results.IM;
d_val = p.Results.d_val;
o_val = p.Results.o_val;
figures = p.Results.figures;
verbose = p.Results.verbose;



if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Segmentation...\n', c(4), c(5), round(c(6))); end
% Whole image dilation. Hopefully this will join some patchy boundaries so that
% objects (tubules) won't leak into each other. The size of the strel can
% probably be tuned depending on the image.
IMd = imdilate(IM, strel('disk', d_val));

% Whole image opening. Should make the boundaries clearer and remove noise coming
% from the interstitial space. strel size should also be adjusted for the
% image (larger strels needed to compensate for greater interstitial space)
IMo = imopen(imcomplement(IMd), strel('disk', o_val));

% Label the objects and overlay on the original image.

% Adaptive thresholding. Segments the image based on local windows of
% maximal intensity
IMbw = imbinarize(IMo, 'adaptive', 'ForegroundPolarity', 'dark');
IMbwl = bwlabel(IMbw);
% Attempt to remove background
    tmp = size(IMbwl);
    tmp2 = [0, IMbwl(1, 1), IMbwl(1, tmp(2)), ...
            IMbwl(tmp(1), 1), IMbwl(tmp(1), tmp(2))];
    tmp3 = unique(tmp2);
    IMbw = ~ismember(IMbwl, tmp3);
% Fill objects with holes in them.
IMbwf = imfill(IMbw, 'holes');
% Convert binary image to connected components with index 1, 2, ...
IMbwl = bwlabel(IMbwf);

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Custom bg removal...\n', c(4), c(5), round(c(6))); end
% Attempt to remove slide as background object and then re-index
IMbwl(IMbwl == IMbwl(1, 1) | ...
    IMbwl == IMbwl(1, size(IMbwl, 2)) | ...
    IMbwl == IMbwl(size(IMbwl, 1), 1) | ...
    IMbwl == IMbwl(size(IMbwl, 1), size(IMbwl, 2))) = 0;
IMbwl = bwlabel(IMbwl);

if figures
    figure
    imshow(IM)
    hold on
    h = imshow(label2rgb(IMbwl, 'jet', 'k'));
    h.AlphaData = 0.3;
end

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Find object boundaries...\n', c(4), c(5), round(c(6))); end
% Compute boundary coordinates of each object, stored in a cell array
% (basically like a list in R). Each cell contains boundary coordinates of
% a single tubule, in order.
IMB = bwboundaries(IMbwl);

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Regionprops...\n', c(4), c(5), round(c(6))); end
% Compute tubule statistics, stored in a struct, which will be used later
% to assign and predict tubule stages.
% R = regionprops(logical(IMbwl), ...
%     'Area', 'Centroid', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', ...
%     'Eccentricity', 'Orientation', 'ConvexArea', 'Circularity', 'FilledArea', ...
%     'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter', ...
%     'PerimeterOld', 'MaxFeretProperties', 'MinFeretProperties');

IMR = regionprops(logical(IMbwl), IM, ...
    'Centroid', 'Area', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', ...
    'Eccentricity', 'Orientation', 'ConvexArea', 'Circularity', 'FilledArea', ...
    'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter', ...
    'PerimeterOld', 'MaxFeretProperties', 'MinFeretProperties', ...
    'MaxIntensity', 'MeanIntensity', 'MinIntensity', 'PixelValues', 'WeightedCentroid');

% Custom params added to regionprops call
for i = 1:length(IMR)
    % R(i).MeanIntensity = mean(int);
    IMR(i).IntDen = IMR(i).MeanIntensity * IMR(i).Area;
    IMR(i).RawIntDen = sum(IMR(i).PixelValues);
end

% Remove PixelValues field (can't save long lists to csv)
IMR = rmfield(IMR, 'PixelValues');