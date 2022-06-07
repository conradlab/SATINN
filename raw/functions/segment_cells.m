function [IMbwl, IMR] = segment_cells(IM, varargin)

%% Segment cells (Hoechst)

% Inputs:
% IM:       Cell image to be segmented.

% Optional inputs, specified as name-value pair arguments:
% method:   Thresholding method for segmentation. Default is 'mean'
%           intensity, 'otsu' might perform better for some of our images.
% hmin:     Value passed to imextendedmin, default is 2.
% figures:  Option to display figures for validation purposes. Default
%           false.
% verbose:  Adds additional information during function runtime. Default
%           false.

p = inputParser;
addRequired(p, 'IM');
addParameter(p, 'method', 'mean', @(x) any(validatestring(x, {'mean', 'custom', 'otsu'})));
addParameter(p, 'hmin', 2, @(x) isnumeric(x));
addParameter(p, 'figures', false, @(x) islogical(x));
addParameter(p, 'verbose', false, @(x) islogical(x));

parse(p, IM, varargin{:});
IM = p.Results.IM;
method = p.Results.method;
hmin = p.Results.hmin;
figures = p.Results.figures;
verbose = p.Results.verbose;



% Alternatively, use a global image threshold (Otsu's method)
if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Segmentation using ', c(4), c(5), round(c(6))); end

if strcmp(method, 'otsu')
    T = graythresh(IM);
    if verbose, fprintf('Otsu: %2.5f\n', T); end
elseif strcmp(method, 'custom')
    T = 0.35;
    if verbose, fprintf('Custom value: %2.5f\n', T); end
else
    T = mean(IM(:))/double(max(IM(:)));
    if verbose, fprintf('Mean: %2.5f\n', T); end
end

xb = imbinarize(IM, T); %x2

% Another alternative, do non-uniform illumination correction before
% segmentation:
% bg = imopen(IM, strel('disk', 30));
% x3 = IM - bg;
% xb = imbinarize(x3);

% You can visualize the differences between these segmented images using
% the imshowpair function, e.g. `figure, imshowpair(x1, x2)`.

% We can't directly apply the watershed function to this binary image
% because we care about ridges within each object (which would indicate
% multiple cells that were merged). Instead, we need to:

% Apply distance transform. Feel free to see what happens to our image
% after each step using `imshow(D)`, etc.
if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Object Splitting...\n', c(4), c(5), round(c(6))); end
D = -bwdist(~xb);

% This next step accounts for oversegmentation. Small minima within objects
% can become a watershed basin, which then results in unnecessary splitting
% of objects into more cells than actually exist. imextendedmin attempts to
% filter out these small minima by looking for "true cell centers"...
mask = imextendedmin(D, hmin);
% (You can tune the H value (the second value you pass to imextendedmin)
% based on whether you want to under- or over-predict cell separation.
% Lower values, e.g. 1 or even 0, will tend to split more uncertain objects
% into multiple cells while higher values, e.g. 3 or 4, will fuse together
% overlapping objects into single cells. In the end, this param will likely
% be image- and even section- or mutant-depdendent.)

% ...while imimposemin modifies the distance transform image (D) to force
% minima to only exist at the centers of those predicted cells.
D2 = imimposemin(D, mask);

% Now we can run watershed on the resulting image...
Dw = watershed(D2);

% (Duplicate binary image)
xbw = xb;
% ... and apply those ridges to original binary image in order to make our
% best guess at splitting overlapping cells.
xbw(Dw == 0) = 0;

% Index cells
IMbwl = bwlabel(xbw);

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Total objects found: %i\n', c(4), c(5), round(c(6)), max(IMbwl(:))); end

if figures
    % Show result
    c = clock;
    fprintf('%02.0f:%02.0f:%02.0f: Drawing figure 1...\n', c(4), c(5), round(c(6)));
    figure
    imshow(IM)
    hold on
    h = imshow(label2rgb(IMbwl, 'jet', 'k'));
    h.AlphaData = 0.3;
    
end

% Skip this if not requesting regionprops
if nargout == 2
    % Calculate regionprops
    if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Regionprops...\n', c(4), c(5), round(c(6))); end
    % Rx = regionprops(logical(xbwl), ...
    %     'Area', 'Centroid', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', ...
    %     'Eccentricity', 'Orientation', 'ConvexArea', 'Circularity', 'FilledArea', ...
    %     'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter', ...
    %     'PerimeterOld', 'MaxFeretProperties', 'MinFeretProperties');
    IMR = regionprops(logical(IMbwl), IM, ...
        'Area', 'Centroid', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', ...
        'Eccentricity', 'Orientation', 'ConvexArea', 'Circularity', 'FilledArea', ...
        'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter', ...
        'PerimeterOld', 'MaxFeretProperties', 'MinFeretProperties', ...
        'MaxIntensity', 'MeanIntensity', 'MinIntensity', 'PixelValues', 'WeightedCentroid');
    
    % Custom params added to regionprops call
    if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Add custom params to regionprops...\n', c(4), c(5), round(c(6))); end
    for i = 1:length(IMR)
        IMR(i).IntDen = IMR(i).MeanIntensity * IMR(i).Area;
        IMR(i).RawIntDen = sum(IMR(i).PixelValues);
    end
    
    % Remove PixelValues field (can't save long lists to csv)
    IMR = rmfield(IMR, 'PixelValues');
    
    % Another validation check fig
    if figures && size(IMR, 1) < 500
        c = clock;
        fprintf('%02.0f:%02.0f:%02.0f: Drawing figure 2...\n', c(4), c(5), round(c(6)));
        b = bwboundaries(IMbwl);
        r2 = struct2table(IMR);
        r2.ID = (1:size(r2, 1))';
        r2(r2.Area <= 5, :) = [];
        col = parula(size(b, 1));
        figure
        imshow(IM)
        hold on
        for i = 1:length(b)
            boundary = b{i};
            plot(boundary(:, 2), boundary(:, 1), 'Color', col(i, :), 'LineWidth', 1);
        end
        for i = 1:size(r2, 1)
            text(r2.Centroid(i, 1), r2.Centroid(i, 2), string(r2.ID(i)), 'FontSize', 10, 'FontName', 'FixedWidth', 'Color', [1 0 0], 'HorizontalAlignment', 'center')
        end
    end
    
end