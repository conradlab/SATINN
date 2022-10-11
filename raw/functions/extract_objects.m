function ims = extract_objects(IM, IMbwl, IMR, annotations, outIMsize, varargin)

% Inputs:
% IM:       Raw image source. Recommend that this be pre-processed as
%           necessary, e.g. using `imtophat` to normalize signal.
% IMbwl:    Segmented image source, e.g. output of segment_cells. This is
%           required if removal of non-specific signal is desired,
%           otherwise, this argument must be empty [], in order to include
%           all non-specific signal within each image.
% IMR:      Regionprops of image source, e.g. output of segment_cells.
%           Specifially, needs the Centroid locations of ROIs.
% annotations: table, with column 1 being cell ID corresponding to the
%           segmented version of image IM, and column 2 being the
%           annotations, e.g. celltypes, as strings. This argument can be
%           empty [] to include all centroids present in IMR.
% outIMsize: specifies the size of the output images, i.e. expected ROI
%           size of each object to be extracted

% Optional inputs, specified as name-value pair arguments
% outdir:   Directory to save extracted objects to (will be created if it
%           does not exist). If this argument is left empty, no images will
%           be saved. The data can be accessed through the output variable
%           `ims` regardless of what is passed to `outdir`.
% figures:  Option to display figures for validation purposes. Default
%           false.
% verbose:  Adds additional information during function runtime. Default
%           false.

% Outputs (optional):
% ims:      Images, stored in a cell array. Additional columns contain
%           annotation data for convenience. Specification of this variable
%           is not required to save images, but if saved images are not
%           desired, the `outdir` argument should be empty.

% Example:
% IM = imread(...)
% IMp = imadjust(imtophat(IM, strel('disk', 100)));
% annotations = readtable(...);
% [IMbwl, ~, IMR] = segment_cells(IMp, ...);
% extract_objects(IM, IMbwl, IMR, annotations, 50, 'outdir', '../outputs/testdir', 'verbose', true);


p = inputParser;
addRequired(p, 'IM');
addRequired(p, 'IMbwl');
addRequired(p, 'IMR', @(x) isstruct(x));
addRequired(p, 'annotations', @(x) istable(x) || iscell(x) || isempty(x));
addRequired(p, 'outIMsize', @(x) isnumeric(x));
addParameter(p, 'outdir', [], @(x) isempty(x) || ischar(x));
addParameter(p, 'verbose', false, @(x) islogical(x));
addParameter(p, 'figures', false, @(x) islogical(x));


parse(p, IM, IMbwl, IMR, annotations, outIMsize, varargin{:});
IM = p.Results.IM;
IMbwl = p.Results.IMbwl;
IMR = p.Results.IMR;
annotations = p.Results.annotations;
outIMsize = p.Results.outIMsize;
outdir = p.Results.outdir;
verbose = p.Results.verbose;
figures = p.Results.figures;


% extract ID'd cells

% image source:
% IM2 = imadjust(imtophat(IM, strel('disk', 100)));

if ~isempty(IMbwl) && (size(IM, 1) ~= size(IMbwl, 1) || size(IM, 2) ~= size(IMbwl, 2))
    error('Mismatch input image sizes (IM and IMbwl).');
end

if isempty(outdir)
    %     outdir = '../outputs/test';
    if verbose, fprintf('No output directory specified.\n'); end
else
    if verbose, fprintf('Saving images to %s..\n', outdir); end
    
    if ~exist(outdir, 'dir')
        if verbose, fprintf('Specified directory (%s) does not exist. Creating...\n', outdir); end
        mkdir(outdir);
    end
    
    flist = dir(fullfile(outdir, '*.tif'));
    %     if strcmp(outdir, '../outputs/test') && ~isempty(flist)
    %         if verbose, fprintf('outdir (likely unspecified) is temporary storage location. All previous files will be removed before creating new ones.\n'); end
    %         for i = 1:length(flist)
    %             delete(fullfile(outdir, flist(i).name));
    %         end
    %     else
    if exist(outdir, 'dir') && ~isempty(flist)
        if verbose, fprintf('Warning: Specified directory already exists and contains files, which may be overwritten.\n'); end
    end
end


% image segmentation (do this outside of function)
% c = clock;
% fprintf('%02.0f:%02.0f:%02.0f: Image segmentation (to customize params, check segment_cells)...\n', c(4), c(5), round(c(6)));
% fprintf('---\n');
% [IMbwl, ~, IMR] = segment_cells(IM2, false);
% fprintf('---\n');
% IMR = struct2table(IMR);

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Image object extraction...\n', c(4), c(5), round(c(6))); end

% checks for presence of segmented image (output of segment_cells), which
% specifies ROIs of each object.
if verbose
    if isempty(IMbwl)
        fprintf('  Segmentation image not provided, extracting full images.\n');
    else
        fprintf('  Segmentation image provided. Extracting within-ROI images only.\n');
    end
end

outIMsize = round(outIMsize/2);

if istable(annotations)
    annotations = table2cell(annotations);
elseif isempty(annotations)
    if verbose, fprintf('  No annotations provided, defaulting to use all objects.\n'); end
    annotations = horzcat(num2cell(1:size(IMR, 1))', repmat({'NA'}, size(IMR, 1), 1)); % dummy annotations
end



% check annotations for duplicates and remove them
[~, idx] = unique(cell2mat(annotations(:, 1)));
if verbose, fprintf('  Found %i duplicate annotations (out of %i) and removed them (if any).\n', ...
        size(annotations, 1)-length(idx), size(annotations, 1)); end
annotations = annotations(idx, :);


ims = cell(size(annotations, 1), 3);

% % old
% errors = [];
% for i = 1:size(annotations, 1)
%     fprintf('  Image extraction at index %i (object %i)...\n', i, annotations{i, 1});
%     src = IM;
%     % remove all signal outside of tubule-specific ROI (comment this to
%     % allow non-specific signal in images)
%     if ~isempty(IMbwl)
%         src(repmat(IMbwl ~= annotations{i, 1}, 1, 1, size(src, 3))) = 0;
%     end
%     try
%         ims{i, 1} = src(round(IMR.Centroid(annotations{i, 1}, 2))-outIMsize+1:round(IMR.Centroid(annotations{i, 1}, 2))+outIMsize, ...
%                         round(IMR.Centroid(annotations{i, 1}, 1))-outIMsize+1:round(IMR.Centroid(annotations{i, 1}, 1))+outIMsize, :);
%     catch
%         try
%             ims{i, 1} = src(round(IMR.Centroid_2(annotations{i, 1}))-outIMsize+1:round(IMR.Centroid_2(annotations{i, 1}))+outIMsize, ...
%                             round(IMR.Centroid_1(annotations{i, 1}))-outIMsize+1:round(IMR.Centroid_1(annotations{i, 1}))+outIMsize, :);
%         catch
%             errors = [errors; annotations{i, 1}];
%             if verbose, fprintf('  Skipping: Subimage exceeded original image bounds at index %i (object %i).\n', i, annotations{i, 1}); end
%         end
%     end
% end

% 210126: add padding to remove possibility of errors

IMp = padarray(IM, [outIMsize, outIMsize]);

if ~isempty(IMbwl)
    IMbwlp = padarray(IMbwl, [outIMsize, outIMsize]);
end

for i = 1:size(annotations, 1)
    if verbose, fprintf('  Image extraction at index %i (object %i)...\n', i, annotations{i, 1}); end
    src =            IMp(round(IMR(annotations{i, 1}).Centroid(2))+1:round(IMR(annotations{i, 1}).Centroid(2))+outIMsize*2, ...
        round(IMR(annotations{i, 1}).Centroid(1))+1:round(IMR(annotations{i, 1}).Centroid(1))+outIMsize*2, :);
    if ~isempty(IMbwl)
        srcbwlp = IMbwlp(round(IMR(annotations{i, 1}).Centroid(2))+1:round(IMR(annotations{i, 1}).Centroid(2))+outIMsize*2, ...
            round(IMR(annotations{i, 1}).Centroid(1))+1:round(IMR(annotations{i, 1}).Centroid(1))+outIMsize*2);
        src(repmat(srcbwlp ~= annotations{i, 1}, 1, 1, size(src, 3))) = 0;
    end
    ims{i, 1} = src;
end


% if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Total subimages skipped: %i (%2.1f%%)...\n', c(4), c(5), round(c(6)), length(errors), length(errors)/size(annotations, 1)); end
ims(:, 2:3) = annotations(:, 1:2);

% remove rows containing empty images and sort
ims(cellfun(@isempty, ims(:, 1)), :) = [];
ims = sortrows(ims, 2); % sort is required to set labels correctly in datastore as they are automatically arranged alphanumerically

if size(IM, 3) == 3
    % write images to files
    if ~isempty(outdir)
        if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Writing images to files...\n', c(4), c(5), round(c(6))); end
        for i = 1:size(ims, 1)
            imwrite(ims{i, 1}, strcat(outdir, sprintf('/%06.0f', ims{i, 2}), '.tif'))
        end
        
        %     % create datastore
        %     if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Create datastore...\n', c(4), c(5), round(c(6))); end
        %     imds = imageDatastore(outdir);
        
        %     % append labels
        %     try
        %         imds.Labels = categorical(ims(:, 3));
        %     catch
        %         try
        %             imds.Labels = categorical(cellfun(@num2str, ims(:, 3)));
        %         catch
        %             fprintf('  Warning: Append labels to datastore failed. Perhaps the data type is incorrect? Skipping...\n');
        %         end
        %     end
        
    end
    
    if figures
        % try to display some images in datastore
        if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Displaying some images for verification...\n', c(4), c(5), round(c(6))); end
        
        figure
        n = min(28, length(annotations));
        perm = randperm(size(ims, 1), n);
        for i = 1:n
            subplot(4, 7, i);
            % imshow(imds.Files{perm(i)});
            imshow(ims{perm(i), 1});
        end
    end
else
    fprintf('Input image is not M-by-N-by-3 (RGB). Skipping writing and displaying samples.\nTo save the data as a Matlab variable anyway, use the output argument.\n');
end

if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Complete.\n', ...
        c(4), c(5), round(c(6))); end

end