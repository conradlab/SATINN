function [YPred, YPredCnf] = cnnValidationFigs(net, dsValidation, varargin)

% Description:
% Custom validation figures to evaluate performance of neural net training.
%
% Inputs:
% net:      Neural network resulting from the trainNetwork function.
% dsValidation: Image datastore containing the validation data (not
%           data used to train the network). This is usually partitioned
%           from the main dataset, e.g. by using splitEachLabel.
%
% Optional inputs:
% show:     Max number of incorrect net predictions to show (default: 8) in
%           a figure (single row). Not recommended to exceed ~12.
% cmap:     Custom colormap that contains at least the number of unique
%           classes, e.g. parula(4). Can be left empty to use a default
%           palette.
%
% Outputs:
% Figures, containing the following: (no variables)
% 1:        A sample of incorrect predictions (mismatch between validation
%           annotation and classification), limited to the value set in
%           'show'. Bar plots in the second row show the distribution of
%           assignment for each class for the corresponding image. Colors
%           highlight the predicted class and expected class.
% 2:        A confusion matrix for all of the validation data. For
%           additional information, see Matlab's `confusionmat` function.
% 3:        A histogram of the highest assignment probability for each
%           image, i.e. the confidence for each classification.
% 4:        A plot of confidence, split by class. Unlike the above plot,
%           this plot attempts to discern confidence between classes (i.e.
%           are certain class assignment less confident than others?)
%
% Example:
% [dsTrain, dsValidation] = splitEachLabel(ds, 0.7, 'randomize');
% layers = [
%   imageInputLayer(50, 50, 3)
%   convolution2dlayer(3, 8, 'padding', 'same', 'stride', 2)
%   ...]
% options = trainingOptions('sdgm', ...)
% net = trainNetwork(dsTrain, layers, options);
% cnnValidationFigs(net, dsValidation);

p = inputParser;
addRequired(p, 'net', @(x) isa(x, 'SeriesNetwork'));
addRequired(p, 'dsValidation', @(x) isa(x, 'matlab.io.datastore.ImageDatastore'));
addParameter(p, 'show', 8, @(x) isnumeric(x));
addParameter(p, 'cmap', [], @(x) isempty(x) || isnumeric(x));

parse(p, net, dsValidation, varargin{:});

net = p.Results.net;
dsValidation = p.Results.dsValidation;
n = p.Results.show;
cmap = p.Results.cmap;



%% classify validation data

[YPred, YPredCnf] = classify(net, dsValidation);
YValidation = dsValidation.Labels;
accuracy = sum(YPred == YValidation) / numel(YValidation);

fprintf('Validation accuracy: %2.2f%%\n', accuracy*100);


%% diagnostics

% show subsample of (in)correct validations
% n = 8 % show this many

    labels = unique(dsValidation.Labels);
    if isempty(cmap)
        cmap = parula(length(labels));
    end
    
tmp = imread(dsValidation.Files{1, 1});
if size(tmp, 3) ~= 1 && size(tmp, 3) ~= 3
    warning('Unable to display non grayscale/RGB images.');
else
    val = min(n, sum(YPred ~= YValidation));

    
    perm = randsample(find(YPred ~= YValidation)', val);
    
    figure
    for i = 1:val
        subplot(2, val, i);
        ff = dsValidation.Files{perm(i)};
        imshow(ff);
        [fp, nm, ~] = fileparts(ff);
        fp2 = strsplit(fp, '\');
        title(sprintf('File: %s %s\nPred: %s\nActual: %s', fp2{end}(1:end-8), nm, YPred(perm(i)), YValidation(perm(i))), 'interpreter', 'none');
        subplot(2, val, val+i);
        b = barh(1, YPredCnf(perm(i), :), 'stacked', 'ShowBaseline', 'off', 'BarWidth', 0.1);
        for j = 1:size(cmap, 1)
            b(j).FaceColor = cmap(j, :);
        end
        tmp = find(~ismember(labels, horzcat(cellstr(YPred(perm(i))), cellstr(YValidation(perm(i))))))';
        for j = tmp
            b(j).FaceColor = 'none';
        end
        % b(labels == cellstr(YPred(perm(i)))).FaceAlpha = 0.5; % dim prediction
        % b(labels == cellstr(YValidation(perm(i)))).FaceAlpha = 1; % true class
        axis off
    end
end

% check confusion matrix

c = confusionmat(YValidation, YPred);
figure
cc = confusionchart(c, unique(dsValidation.Labels));
cc.RowSummary = 'row-normalized';
cc.ColumnSummary = 'column-normalized';
title('Confusion matrix');

% check confusion matrix for entries with confidence > 80%

c = confusionmat(YValidation(max(YPredCnf, [], 2) > 0.8), YPred(max(YPredCnf, [], 2) > 0.8));
figure
cc = confusionchart(c, unique(dsValidation.Labels));
cc.RowSummary = 'row-normalized';
cc.ColumnSummary = 'column-normalized';
title(sprintf('Confusion matrix with >80%% confidence\n (Kept %2.1f%% data)', sum(max(YPredCnf, [], 2) > 0.8) * 100 / size(YPredCnf, 1)));

% plot cell assignment confidence

[v, id] = max(YPredCnf, [], 2);

b = 0:0.05:1;
h1 = histcounts(v(dsValidation.Labels == YPred), b);
h2 = histcounts(v(dsValidation.Labels ~= YPred), b);

figure
b1 = bar(h1+h2, 'FaceColor', [0.85 0.325 0.098], 'EdgeColor', 'none');
b1.XData = b(1:end-1) + mean([b(1), b(2)]);
hold on
b2 = bar(h1, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');
b2.XData = b1.XData;
legend({'Incorrect', 'Correct'}, 'location', 'northwest');
xlabel('Confidence');
ylabel('Count');
title('Assignment confidence pooled');

% split based on cell types
figure
hold on
for i = 1:size(YPredCnf, 2)
    [n, e] = histcounts(v(id == i), 0:0.05:1);
    plot(0.5 * (e(1:end-1) + e(2:end)), smooth(n), 'Color', cmap(i, :), 'LineWidth', 2)
end
set(gca, 'YScale', 'log');
axis([0 1 1 inf]);
% create custom legend
classes = unique(YPred);
dummy = zeros(length(classes), 1);
for i = 1:length(classes)
    dummy(i) = plot(NaN, NaN, 'Marker', 'o', 'MarkerEdgeColor', cmap(i, :), 'MarkerFaceColor', cmap(i, :), 'LineStyle', 'none');
end
legend(dummy, classes, 'TextColor', 'k', 'location', 'northwest');
xlabel('Confidence');
ylabel('Count');
title('Assignment confidence by class');
% legend('boxoff');
