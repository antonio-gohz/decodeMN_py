function cohResults = MN_coherence(spikeTrains,varargin)
% improvements: splits the entire amount of decomposed MNs into 2 groups
% and computs coherence 
% based on previous analysis (testCoherenceAcrossMUsinCST.m) the minimum 
% number of MNs to have significant peaks in the coherence estimation is 4
% in each train i.e. 8 decomposed MUs 
% todo: take data instead of path (easier to debug)
% explore 2 ways: 
% 1. fixing number of groups (i.e., splitting into groups of 2)
% 2. fixing number of MUs in each CST (4 minimum) to compare across subjects 
% for either we want random permutations forming 2 groups
% the max number of permutations for a minum of 8 decomposed MUs is equal
% to 8! (40320) from n!(n-k)! where n=k so n!*0! = n!
% we calculate only a max of 100
% close all

% Input erros 
if nargin < 1
   error('Provide spike trains'); 
end

% Set default values
fs=2048;
W=[]; %num of MUs per group
overlap=0;
LW=fs/2; % length pwelch
flagFischer = 0;
% default plot parameters
colors = {'#4477AA'}; 
lineSpecs = '-';
plotFlag=0;

% Process optional input parameters
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i+1};
    switch param
        case 'numMUsperGroup'
            W = value;
        case 'LW'
            LW = value;
        case 'FlagFischerZtrans'
            flagFischer = value;
        case 'overlap'
            overlap = value;
        case 'fs'
            fs = value;
        case 'Color'
            colors = value;
        case 'lineSpecs'
            lineSpecs = value;
        case 'PlotFlag'
            plotFlag = value;
        otherwise
            error('Invalid optional input parameter: %s', param)
    end
end
if size(spikeTrains,1)>size(spikeTrains,2)
    spikeTrains=spikeTrains';
end
if isempty(overlap)
   overlap0=50;  % equivalent to 50%
else
    overlap0=overlap;
end

if ischar(colors)
colors = hex2rgb(colors);
end
L = round((length(spikeTrains)-LW)/(LW*(1-overlap0/100)))+1;
%L=round(length(firing)/(LW))*2; % I think *2 is wrong
COF = 1- (1-0.95)^(1/(L-1));
if flagFischer
COF = sqrt(2*L)*atanh(sqrt(COF)); %fisher's z transformation to test significancy of coherence values.
end
N = size(spikeTrains,1);
if isempty(W)
    W=floor(N/2);
end
pairedCombs = subgroupPairs(N,W); % computes all non repeated pairs
NPAIRS=length(pairedCombs);

[COH,F] = mscohere(detrend(sum(spikeTrains(pairedCombs{1}(1,:),:),1),0),...
    detrend(sum(spikeTrains(pairedCombs{1}(2,:),:),1),0),hanning(round(LW)),overlap,10*fs,fs);
    if flagFischer 
COH = sqrt(2*L)*atanh(COH);
    end
COH = [COH,zeros(length(COH),NPAIRS)]; % creates a matrix for the rest of the COHs
for p = 2:NPAIRS
    h= waitbar(p/NPAIRS);
    [COH(:,p)] = mscohere(detrend(sum(spikeTrains(pairedCombs{p}(1,:),:),1),0),...
        detrend(sum(spikeTrains(pairedCombs{p}(2,:),:),1),0),hanning(round(LW)),overlap,10*fs,fs);
    %keyboard
    if flagFischer 
    COH(:,p) = sqrt(2*L)*atanh(COH(:,p));
    end
end
close(h);
coh_mean=nanmean(COH,2);
coh_std = nanstd(COH,1,2);
if plotFlag
hold on,plot(F,coh_mean,lineSpecs,'Linewidth',1.5,'Color',colors); % t=1 = Wout stim (blue)
jbfill(F',(coh_mean-coh_std)' ,(coh_mean+coh_std)',colors,'none',1,0.3,'off');
hold on,plot(F,ones(1,length(F))*COF,':','Linewidth',1,'Color',colors);
xlim([0,60])
end
cohResults.confThres = COF;
cohResults.coh_mean = coh_mean;
cohResults.pairedCombs = pairedCombs;
cohResults.numMUsperGroup=W;
cohResults.COHstd = coh_std;
cohResults.LW=LW;
cohResults.overlap=overlap0;
cohResults.f= F;

cohResults.peaks=[max(coh_mean(F<5)), max(coh_mean(F>5 & F<15)),max(coh_mean(F>15))]; % delta, alpha, beta
cohResults.areas=[sum(coh_mean(F<5)), sum(coh_mean(F>5 & F<15)),sum(coh_mean(F>15))]; % delta, alpha, beta






