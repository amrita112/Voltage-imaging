function forAmrita
%temporal filtering of signals with a template

resp = exp(-abs(linspace(-8,8,21))).* sin(linspace(pi,-pi, 21)+0.3);
resp = resp./sqrt(mean(resp.^2));

L = 100000; %length of simulation
spikesI = rand(1,L)>0.995; %initial description of spikes; we'll jitter these at superresolution for more realism
spikesI(1) = false;
T = find(spikesI); A = rand(size(T));
spikes = double(spikesI);
spikes(T-1) = A; spikes(T) = 1-A; %jitter spikes

%generate spectrally varying noise
noiseFFT = fft(randn(1,L));
noiseFFT = noiseFFT.*(smooth(rand(1,L)-0.5, 500).^2)';
B = noiseFFT;
B(1001:end-999) = 0;
noiseFFT([1:1000 end-998:end]) = 0;
noise = real(ifft(noiseFFT));
b = real(ifft(B));
b = b./sqrt(mean(b.^2));
noise = noise./sqrt(mean(noise.^2));

%simulate data
amp = 1.2;  %signal-to-noise ratio
GTdata =  conv(amp*double(spikes),resp, 'same');
data=  conv(amp*double(spikes),resp, 'same') + noise + b;

figure('name', 'raw simulated data')
plot(data)

%detrend data
DATA = fft(data, 2*L-1);
nFreqNull = round(L/100);
DATA([1:nFreqNull end-nFreqNull+2:end]) = 0;
data = ifft(DATA); 
data = data(1:L);


dataS = data';
%dataS = smooth(data,3);%smooth data for initial detection?

%threshold
figure, hist(findpeaks(dataS), 1000)
SIGMA = std(dataS(dataS<prctile(data,99)));
thresh = min(prctile(dataS, 99.99), 3.5*SIGMA);

%thresh = setThresh(findpeaks(dataS), 0.999);
[~,locs] = findpeaks(dataS, 'MinPeakHeight',thresh, 'MinPeakDistance', 4);
figure('name', 'initial spike detection'), plot(dataS), hold on, scatter(locs, thresh*ones(size(locs)))

%peak-triggered average
window = -10:10;
locs = locs(locs>(-window(1)+1) & locs<(length(dataS)-window(end)));
PTD = data(locs+repmat([-10:10], size(locs,1),1));
PTA = mean(PTD,1);%peak-triggered average
%PTA = PTA-mean(PTA);
%PTA = PTA.*(2*min(hanning(length(PTA)),0.5))'; %window it so it decays to 0, important for fft later

%plot how well we're doing in simulation
locsGT = find(spikesI)';
locsGT = locsGT(locsGT>(-window(1)+1) & locsGT<(length(dataS)-window(end)));
PTDgt = data(locsGT+repmat(window, size(locsGT,1),1));
figure('name', 'Spike triggered averages'), plot(PTA), hold on, plot(mean(PTDgt,1)), plot(resp*amp)
legend({'Peak Triggered Average', 'Peak Triggered Average- Ground Truth', 'True Spike'})

%wiener filter
% Nfft = 2.^nextpow2(length(data))-1;
% FFTdata = fft(data, Nfft);

guessdata = zeros(size(data));
guessdata(locs) = 1;
guessdata = conv(guessdata,PTA, 'same');
noisedata = data-guessdata;
% FFTnoise= fft(noisedata, Nfft);

% guessNspikes = 10*length(locs); %to do: use histogram to estimate what fraction of spikes we initially detected
% FFTexpected = nan(500,Nfft);
% for iter = 1:500
%     SS = double(rand(size(data))<(guessNspikes/numel(data))); SS(1) = 0;
%     TT = find(SS & ~[0 SS(1:end-1)]); A = rand(size(TT));
%     SS(TT-1) = A; SS(TT) = 1-A; %jitter spikes
%     idealdata = conv(SS,PTA, 'same');
%     FFTexpected(iter,:) = real(fft(idealdata)).^2;
% end
% S2 = mean(FFTexpected,1);
% datafilt = wienerFilterKP(S2,FFTdata,FFTnoise);
%pxx = pwelch(x,window,noverlap,nfft)

%matched filter
%prewhiten
datafilt = whitenedMatchedFilter(data, locs, window);

%datafilt2 = wienerFilter(guessdata,data,noisedata);

%spikes detected after filter
SIGMA = std(datafilt(datafilt<prctile(data,99)));
thresh2 = min(prctile(datafilt, 99.99), 3.5*SIGMA);
[~,locs2] = findpeaks(datafilt, 'MinPeakHeight',thresh2, 'MinPeakDistance', 4);

P1 = performance(locs, locsGT); P2 = performance(locs2, locsGT);


figure('name', 'performance'); plot(P1); hold on, plot(P2)
set(gca, 'xtick', [1 2 3], 'xticklabel', {'true Positive', 'false Positive', 'misses'})
legend({'Pre', 'Post'});
figure('name', 'histograms original vs filt'), 
subplot(2,1,1); hist(findpeaks(data),1000); hold on, plot([thresh thresh], get(gca, 'ylim'))
subplot(2,1,2);  hist(findpeaks(datafilt), 1000); hold on, plot([thresh2 thresh2], get(gca, 'ylim'))

% normalize and display result
datafilt= datafilt./norm(datafilt);
data = data./norm(data);
GTdata = GTdata./norm(GTdata);
figure('name', 'traces original vs filt'), 
ax1 = subplot(2,1,1); plot(data); hold on, scatter(locs, max(data)*ones(size(locs))); scatter(find(spikesI), max(data)*1.1*ones(1,sum(spikesI)), '*');
ax2 = subplot(2,1,2); plot(datafilt); hold on, scatter(locs, max(datafilt)*ones(size(locs))); scatter(find(spikesI), max(datafilt)*1.1*ones(1,sum(spikesI)),'*'); scatter(locs2, max(datafilt)*0.9*ones(size(locs2)))
legend({'trace', 'original detections', 'GT', 'postfilter'})
linkaxes([ax1 ax2]);

newguessdata = zeros(size(data));
newguessdata(locs2) = 1;
newguessdata = conv(newguessdata,PTA, 'same');

disp('Correlation of data to Ground truth:')
corr(GTdata(:), data(:))
disp('Correlation of filtered data to Ground truth:')
corr(GTdata(:), datafilt(:))
disp('Correlation of spike times to Ground truth:')
corr(GTdata(:), newguessdata(:))
% disp('Correlation of filtered data to Ground truth:')
% corr(GTdata', datafilt1')
end

function P = performance(locs, locsGT)
%spikes detected originally
    distMTX = locs(:)-locsGT(:)';
    truePos = sum(min(abs(distMTX),[],2)<=2);
    falsePos = size(distMTX,1)-truePos;
    misses = size(distMTX,2)-truePos;
    P = [truePos falsePos misses];
end
% %detrend
% N = fft(noise);
% LF = 100;
% N([1:LF end-LF+1:end]) = 0;
% noiseDT = real(ifft(N));