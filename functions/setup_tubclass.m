function [im, imabwl, imaR] = setup_tubclass(imh, ima, imv)

% import data

% normalize signal intensity
imhp = imadjust(imtophat(imh, strel('disk', 100)));
imap = imadjust(imtophat(ima, strel('disk', 100)));
imvp = imadjust(imtophat(imv, strel('disk', 100)));
imvp = imadjust(imvp, [0.1, 1]); % specifically for acrv1, filter out the bottom 10% of signal as bg post-normalization (empirically determined)

if nargout > 1
% segment acta2 channel
imabwl = segment_tubules(imap);

% tubule ROI dilation (compensate for interstitial space removal)
imabwld = imdilate(imabwl, strel('disk', 30));

% tubule ROI statistics
imaR = regionprops(imabwld, {'Area', 'Centroid', 'Circularity'});

% automatically filter out unlikely tubules (over/under-segmentation, etc.
% based on statistical criteria; parameters may differ depending on the
% image used)
filter = [imaR.Circularity] > 0.5 & [imaR.Area] > 0.5*10^6 & [imaR.Area] < 3*10^6;
imabwl(~ismember(imabwl, find(filter))) = 0; % remove objects that don't meet criteria
% reindex and reestablish regionprops (optional, but useful to prevent
% extracting a lot of noise/non-tubules, which would have to be filtered
% out anyway during manual annotation)
imabwld = imdilate(bwlabel(imabwl), strel('disk', 30));
imaR = regionprops(imabwld, {'Area', 'Centroid', 'Circularity'});
end

% assemble RGB image
im = cat(3, imvp, imhp, imap);
