function [signalRef, decompParametersRef] = ...
    editMUFilters(signalRef, decompParametersRef,EMGs,fs,manual_edit,plotFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization loop of the MU filter to minimize the Silhoutte

% Input: 
%   w = initial weigths
%   X = whitened signal
%   SIL = Silhoutte
%   fsamp = sampling frequency


% Output:
%   wlast = new weigths (MU filter)
%   spikeslast = discharge times of the motor unit
%   SILlast = coefficient of varation of the inter spike intervals 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%EMGs = EMGFilt(ids_decomp1(1):ids_decomp1(2),:)';
%manual_edit=false; initial commit
t= linspace(0,length(EMGs)/fs,length(EMGs))';

X = extend(EMGs(decompParametersRef.EMGmask,:), decompParametersRef.extensionfactor);
X = X(:,1:size(EMGs, 2));
%X = decompParametersRef.whitMat*X;
%X = decompParametersRef.normtanh.*tanh(X./decompParametersRef.normtanh);
if nargin<6
    plotFlag = false;
end

if manual_edit % format plots 
    formats_init = {'.','-','o'};
        plotFlag = true;
else
    formats_init = {'xr','-r','or'};
end



idMU=1;
while idMU<=size(decompParametersRef.MUFilters,2)
    % init params
    wini = decompParametersRef.MUFilters(:,idMU);
    w =decompParametersRef.MUFilters(:,idMU);
    spikes = logical(signalRef.spikeTrains(idMU,:));
    ipt = signalRef.Pulsetrain(idMU,:);
    drs= signalRef.dischargeRates(idMU,:);
    normIPT=decompParametersRef.normIPT(idMU);
    iReSIGt=decompParametersRef.iReSIGt;
    if plotFlag
    % init plots
    figure('units','normalized','Position',[0 0 0.5 1]);
    tl= tiledlayout(4,1);
    ax(1) = nexttile;
    plot(t,drs,formats_init{1},'MarkerSize',4)
    ylabel('Discharge rates')
    ax(2) = nexttile;
    plot(t,spikes,formats_init{2})
    ylabel('Spike trains' )

    ax(3) = nexttile;
    plot(t,ipt,formats_init{2},t(spikes),ipt(spikes),formats_init{3})
    ylabel('Pulse trains')

    linkaxes(ax(1:3),'x')
    xlim([min(t(spikes))-0.5,max(t(spikes))+0.5]);
    sgtitle("MU: "+idMU + ". Left clic: choose spike to eliminate // Right clic: next MU // Middle clic: Ends all")
    ax(4) = nexttile;
    ylabel('MU filter')
    end
    if manual_edit
        hfil = plot(1:length(wini),wini,1:length(wini),w);
        hold(ax(1:3),'on')

        while(1)
            [Imx, Imy , mb]=ginput(1);
            if mb==1 % left clic -> eliminates points
                % find closest spike to users choice (in any subplot as it is
                % referenced to x axis
                [~, minid] = min(abs(Imx-t(spikes)));
                ind = find(spikes,minid,'first');
                ind= ind(end);

                %updatePlot(ax,t,drs,ind,w)
                % plot selected markers
                plot(ax(1),t(ind),drs(ind),'.r','MarkerSize',10)
                plot(ax(2),t(ind-1:ind),spikes(ind-1:ind),'-r')
                plot(ax(3),t(ind-2:ind+2),ipt(ind-2:ind+2),'-r')
                xlim([min(t(spikes))-0.5,max(t(spikes))+0.5]);

                spikes(ind-1:ind) =  0;
                drs(ind) = nan;

                w = mean(X(:,spikes), 2);
                [ipt, distTimes,sil, normIPT,centroids] = getspikes(w, X, fs, iReSIGt);
                spikes=false(size(spikes));
                spikes(distTimes) = true;
                plot(ax(3),t(ind-2:ind+2),ipt(ind-2:ind+2),'-b')

                set(hfil,{'XData'},{(1:length(wini))',(1:length(wini))'}',{'YData'},{wini,w}')
                %plot(ax(4),1:length(wini),wini,1:length(wini),w)
            end
            if mb==3 % next
                answer = questdlg('Save changes?', 'Save changes?',...
                    'Yes','No (do it again)','Yes');
                % Handle response
                switch answer
                    case 'Yes'
                        decompParametersRef.MUFilters(:,idMU)=w;
                        signalRef.spikeTrains(idMU,:)= spikes;
                        signalRef.Pulsetrain(idMU,:)= ipt;
                        signalRef.dischargeRates(idMU,:) =drs;
                        signalRef.SIL(idMU)= sil;
                        COV = std(drs,'omitnan' )/mean(drs,'omitnan' );
                        signalRef.Dischargetimes{idMU} = distTimes;
                        signalRef.COV(idMU) = COV;
                        decompParametersRef.normIPT(idMU) = normIPT;
                        decompParametersRef.centroid(idMU,:) = centroids;
                        %close gcf
                        idMU = idMU+1;
                end
                break;
            end
            if mb==2 % Middle mouse button breaks all
                %close gcf
                idMU = size(decompParametersRef.MUFilters,2)+1;
                break;
            end
        end
        %% automatic editing based on outliers
    else
        if plotFlag
        plot(ax(4),1:length(w),w,formats_init{2});
        hold(ax,"on")
        end
        sp_ids= find(spikes);
        out_ids=find(isoutlier(drs));
        removed_ind =[];
        gradd=1;
        while length(out_ids)>1 && gradd>0.1
            for i=1:length(out_ids)
                [~,ids]= mink(abs(sp_ids-out_ids(i)),3);
                [~,idmin_ica]=min(ipt(sp_ids(ids)));
                ind=sp_ids(ids(idmin_ica));
                removed_ind= [removed_ind,ind];
                spikes(ind-1:ind) =  0;
            end
            drs = nan(size(spikes));
            drs(spikes) = [nan,fs./diff(find(spikes))];
            prev_len = length(out_ids);
            out_ids=find(isoutlier(drs));
            gradd= abs(length(out_ids)-prev_len)./prev_len;
        end
        w = mean(X(:,spikes), 2);
        [ipt, distTimes,sil, normIPT,centroids] = getspikes(w, X, fs, iReSIGt);
        spikes=false(size(spikes));
        spikes(distTimes) = true;
        if plotFlag
        plot(ax(1),t,drs,'.b','MarkerSize',12)
        plot(ax(2),t,spikes,'-b')
        plot(ax(3),t,ipt,'-b',t(spikes),ipt(spikes),'ob')
        plot(ax(4),1:length(w),w,'-b')
        legend({'Before editing','After editing'})
        end
        %update
        COV = std(drs,'omitnan' )/mean(drs,'omitnan' ); 
        signalRef.Dischargetimes{idMU} = distTimes;
        signalRef.COV(idMU) = COV;
        signalRef.spikeTrains(idMU,:)= spikes;
        signalRef.Pulsetrain(idMU,:)= ipt;
        signalRef.dischargeRates(idMU,:) =drs;
        signalRef.SIL(idMU)= sil;
        decompParametersRef.MUFilters(:,idMU)=w;
        decompParametersRef.normIPT(idMU) = normIPT;
        decompParametersRef.centroid(idMU,:) = centroids;
        idMU = idMU+1;
    end

end