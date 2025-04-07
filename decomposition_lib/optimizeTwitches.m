function [apsLast, tcsLast, neuralAct] = ...
    optimizeTwitches(spikeTrains, torque, aps, tcs,showPlots)

 % compute neuralAct
    [neuralAct, ~] = muActivation(spikeTrains,[],...
        [], aps, (tcs)./1000, [], 'Fuglevand',...
        [], [], 0, [],[], [] ,1,[], ...
        [],[],[], [],0);

    [corr,lags] = xcorr(detrend(sum(neuralAct,2)),detrend(torque),'normalized');

    [r,idmax] = max(corr);
    l = lags(idmax);
mean0=mean(aps);
cov0 = std(aps)/mean0;
meantcs0 = mean(tcs);
k = 1;
strategy = 1; % increase sparseness
changeStrategy = 0;
conditionAps = true;
conditionTcs = true;
if showPlots
        %axes properties
    colors =[ 0.2 0.2 0.2;0.2667    0.4667    0.6667; 0.9333    0.4000    0.4667;
              0.2667    0.4667    0.6667; 0.9333    0.4000    0.4667;
              0.2667    0.4667    0.6667; 0.9333    0.4000    0.4667];
    ylimits = [0,1;0,1;0.45*min(tcs),1.1*max(tcs)];
    xlimits = [1,2048*10;0,length(aps);0,length(tcs)];
    ylabels = {'Neural Activation (a.u.)','Peak amplitude(a.u.)','Contraction time (ms)'};
    xlabels = {'Samples','MU (#)','MU (#)'};
    lineStyles = {':','-','-','-','-','-','-'};
    markerStyles = {'none','none','none','o','o','o','o'};
    axesLineConfig = [1,1;    % subplot 1 -signal 1 (IPT)
        1,2; 1,3;      % subplot 1 -signal 2 (IPT peaks)
        2,4; 2,5;
        3,6;3,7];      % % subplot 2 -signal 3 (MU filter)
    [h,ax] = livePlots.preparePlots(1,'linewidth',1,'colors',colors,...
        'ylimits',ylimits ,'xlimits',xlimits,'xlabels',xlabels,'ylabels',ylabels,...
        'lineStyles',lineStyles,'markerStyles',markerStyles,...
        'axesLineConfig',axesLineConfig);
    [h,ax] = livePlots.updatePlot([{1:length(torque),torque'/max(torque)};{1:length(neuralAct),sum(neuralAct,2)'/max(sum(neuralAct,2))};{1:length(neuralAct),sum(neuralAct,2)'/max(sum(neuralAct,2))}; ...
        {1:length(aps),aps'};{1:length(aps),aps'}; ...
        {1:length(tcs),tcs'};{1:length(tcs),tcs'}],...
        h,ax);
    title(ax(1),"Optmization twitches It: "+k);
end

while conditionAps && conditionTcs && k<100
    rlast=r;
    llast = l;
    strategyLast= strategy;
    apsLast=aps;
    tcsLast =tcs;
    % update aps
    aps_norm = aps;%4*(aps-mean0)/(max(aps-mean0))-2;
    if strategy ==1
        aps= aps_norm.*aps_norm;%(aps_norm).*abs(aps_norm); % increases sparseness
    else
        aps=0.125*sqrt(abs(aps_norm)).*sign(aps_norm); % decreases sparseness
    end
    % aps=(aps+2)*max(aps-mean0)/4+mean0; % returns to same scale
    aps = tanh(aps*mean0/mean(aps));
    %[aps,aps2,aps3]

    % update tcs
    tcs = tcsLast-0.08*mean(tcsLast);

    % compute neuralAct
    [neuralAct, ~] = muActivation(spikeTrains,[],...
        [], aps, (tcs)./1000, [], 'Fuglevand',...
        [], [], 0, [],[], [] ,1,[], ...
        [],[],[], [],0);

    [corr,lags] = xcorr(detrend(sum(neuralAct,2)/max(sum(neuralAct,2))),detrend(torque/max(torque)),'normalized');

    [r,idmax] = max(corr);
    l = lags(idmax);
    if (std(aps)/mean(aps)>2 *cov0 || std(aps)/mean(aps)<0.5*cov0) && ~changeStrategy
        changeStrategy = 1;
    else 
        changeStrategy = 0;        
    end
    if ~changeStrategy %r > rlast  % continue strategy
        strategy = strategyLast;
    else   % change strategy
        strategy = setdiff( [1,2],strategyLast);
    end

    conditionAps = abs(r - rlast)>0.0001 & ~any(aps<=0);
    conditionTcs = abs(l - llast)/llast>0.0001 & ~any(mean(tcs)<=0.4*meantcs0);
    k=k+1;
    if showPlots
    [h,ax] = livePlots.updatePlot([{[],[]};{[],[]};{1:length(neuralAct),sum(neuralAct,2)'/max(sum(neuralAct,2))}; ...
        {[],[]};{1:length(aps),aps'}; ...
        {[],[]};{1:length(tcs),tcs'}],...
        h,ax);
    title(ax(1),"Optmization twitches It: "+k+ " Corr = " +r + " delay = "+100*l/2048+" ms");
    title(ax(2),"Strategy "+strategy)
    pause(0.1)
    end
end
