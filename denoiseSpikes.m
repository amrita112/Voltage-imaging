function [datafilt spikeTimes guessData] = denoiseSpikes(data, windowLength)

doSmooth = true;
if doSmooth
    dataS = ((data + data([2:end 1]))/2)';%smooth data for initial detection
else
    dataS = data';
end

%threshold
SIGMA = std(dataS(dataS<prctile(data,99)));
thresh = min(prctile(dataS, 99.99), 3.5*SIGMA);

figure('name', 'Amplitude histogram'), hist(findpeaks(dataS), 1000), hold on, plot([thresh thresh], get(gca, 'ylim'), 'r')
disp(['Threshold set to ' num2str(thresh)]); drawnow;

[~,locs] = findpeaks(dataS, 'MinPeakHeight',thresh, 'MinPeakDistance', 4);

if doSmooth
    %fix aliasing from initial smoothing
   for ix = 1:length(locs)
       [~, maxind] = max(data(locs(ix)+[0:1]));
       locs(ix) = locs(ix) + maxind-1;
   end
end
    
%peak-triggered average
window = -windowLength:windowLength;
locs = locs(locs>(-window(1)+1) & locs<(length(dataS)-window(end)));
PTD = data(locs+repmat(window, size(locs,1),1));
PTA = mean(PTD,1); %peak-triggered average
PTA = PTA-mean(PTA([1 end]));
%PTA = PTA.*(2*min(hanning(length(PTA)),0.5))'; %window it so it decays to 0, important for fft
figure('name', 'Peak-triggered average'), plot(PTD', 'color', [0.5 0.5 0.5]), hold on, plot(PTA, 'k', 'linewidth', 2)

%matched filter
datafilt = whitenedMatchedFilter(data,locs, window);

%wiener filter
% idealdata = zeros(size(data));
% idealdata(locs) = 1;
% idealdata = conv(idealdata,PTA, 'same');
% noisedata = data-idealdata;

% idealdata = zeros(size(data));
% idealdata(locs) = 1;
% idealdata = conv(idealdata,PTA, 'same');
% noisedata = data-idealdata;

% idealdata = zeros(size(data));
% idealdata(randsample(length(idealdata), 4*length(locs))) = 1;
% idealdata = conv(idealdata,PTA, 'same');

%wiener filter
%datafilt2 = wienerFilter(idealdata,data,noisedata);

%spikes detected after filter
SIGMA = std(datafilt(datafilt<prctile(data,99)));
thresh2 = min(prctile(datafilt, 99.99), 3.5*SIGMA);
[~,spikeTimes] = findpeaks(datafilt, 'MinPeakHeight',thresh2, 'MinPeakDistance', 4);

guessData = zeros(size(data));
guessData(locs) = 1;
guessData = conv(guessData,PTA, 'same');

%wiener filter shrinks the data;
%rescale so that the mean value at the peaks is same is in the input
datafilt = datafilt.*(mean(data(locs))./mean(datafilt(locs)));

figure('name', 'Data before and after filtering')
ax1 = subplot(2,1,1);
plot(data), hold on, scatter(locs, max(datafilt)*ones(size(locs)), 'markeredgecolor', 'g')
ax2 = subplot(2,1,2);
plot(datafilt), hold on, scatter(locs, max(datafilt)*ones(size(locs)), 'markeredgecolor', 'g')
linkaxes([ax1 ax2]);
end