function [im, imabwd, imaR] = setup_tubclass(imh, ima, imv)

% import data

% normalize signal intensity
imhp = imadjust(imtophat(imh, strel('disk', 100)));
imap = imadjust(imtophat(ima, strel('disk', 100)));
imvp = imadjust(imtophat(imv, strel('disk', 100)));
imvp = imadjust(imvp, [0.1, 1]); % specifically for acrv1, filter out the bottom 10% of signal as bg post-normalization (empirically determined)

if nargout > 1
    % segment acta2 channel
    % imabwl = segment_tubules(imap);
    
    % if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Segmentation...\n', c(4), c(5), round(c(6))); end
    % Whole image dilation. Hopefully this will join some patchy boundaries so that
    % objects (tubules) won't leak into each other. The size of the strel can
    % probably be tuned depending on the image.
    IMd = imdilate(imap, strel('disk', 25));
    
    % Whole image opening. Should make the boundaries clearer and remove noise coming
    % from the interstitial space. strel size should also be adjusted for the
    % image (larger strels needed to compensate for greater interstitial space)
    IMo = imopen(imcomplement(IMd), strel('disk', 100));
    
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
    IMbwl = imfill(IMbw, 'holes');
    % Convert binary image to connected components with index 1, 2, ...
    % IMbwl = bwlabel(IMbwf);
    
    % if verbose, c = clock; fprintf('%02.0f:%02.0f:%02.0f: Custom bg removal...\n', c(4), c(5), round(c(6))); end
    % Attempt to remove slide as background object and then re-index
    IMbwl(IMbwl == IMbwl(1, 1) | ...
        IMbwl == IMbwl(1, size(IMbwl, 2)) | ...
        IMbwl == IMbwl(size(IMbwl, 1), 1) | ...
        IMbwl == IMbwl(size(IMbwl, 1), size(IMbwl, 2))) = 0;
    % IMbwl = bwlabel(IMbwl);
    
    % tubule ROI dilation (compensate for interstitial space removal)
    imabwd = bwmorph(IMbwl, 'thicken', 40);
    
    % tubule ROI statistics
    imaR = regionprops(imabwd, {'Area', 'Circularity'});
    
    % automatically filter out unlikely tubules (over/under-segmentation, etc.
    % based on statistical criteria; parameters may differ depending on the
    % image used)
    filter = [imaR.Circularity] > 0.5 & [imaR.Area] > 0.5*10^6 & [imaR.Area] < 3*10^6;
    imabwd = bwlabel(imabwd);
    imabwd(~ismember(imabwd, find(filter))) = 0; % remove objects that don't meet criteria
    % reindex and reestablish regionprops (optional, but useful to prevent
    % extracting a lot of noise/non-tubules, which would have to be filtered
    % out anyway during manual annotation)
    
    imabwd = logical(imabwd);
    imaR = regionprops(imabwd, {'Area', 'BoundingBox', 'Centroid', ...
        'Circularity', 'Eccentricity', 'MajorAxisLength', ...
        'MinorAxisLength', 'Orientation', 'Perimeter'});
    
end

% assemble RGB image
im = cat(3, imvp, imhp, imap);