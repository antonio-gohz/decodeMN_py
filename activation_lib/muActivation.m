function [neuralActivation, normSat] = muActivation(spikeTrains, dischargeRate, t, Ap, tc, thr, twitchModel, satLevel, tetParam, tetanusModel, c1,c2, shapeFactor_ ,activationScale_, force,...
    muscle,subject,position, activation,normalization)

    neuralActivation = zeros(size(spikeTrains));
    neuralMatrix = zeros(max(sum(spikeTrains)),size(spikeTrains,2));
    n_spike = 0;
    firstSpike=0;    
    t0s = NaN(max(sum(spikeTrains)),size(spikeTrains,2));
    %t0s = zeros(max(sum(excitation)),size(excitation,2));
    fs=2048;
    T = 1/fs;
    normSat = 1/size(spikeTrains,2).* satLevel;

   % thr = tc.*2.3;
    %% Tetanus parameters
    if tetanusModel
        %tetParam=100; % the bigger ther smoother saturation
        normSat = 1/size(spikeTrains,2).* satLevel;
    end
    %% Raikova parameters
    if strcmp(twitchModel,'Raikova')
        k=log(2)./(thr-tc-tc.*log(thr./tc));
        m =k.*tc;
        p= Ap.*exp(-k.*tc.*(log(tc)-1));
    end
    %% CEINMS variables 
    beta1_ = c1 + c2;
    beta2_ = c1 .* c2;
    alpha_ = 1 + beta1_ + beta2_;
    neuralActivation_=neuralActivation;

    %% test normalization method
%         [tStamps] = extractTimeStamps(muscle,subject,position, activation);
%     numMUs = [];
%     fs = 2048;
%     try
%     for i=1:size(tStamps,1)
%         numMUs = [numMUs,sum(sum(excitation(floor(fs*tStamps(i,1)):floor(fs*tStamps(i,2)),:))>0)];
%     end
%     catch 
%     numMUs=1;
%     end
%     % This is a small band aid very large differences of the number of MUs
%     %make the normalization amplityde very different e.g. 1/5 vs 1/30 MUs
%     % this happens a lot in the soleau theregore limit the numMUs to a
%     % minimum of 10 MUs
%     numMUs(numMUs<10)=10;
%     numMUs(numMUs>30)=30;
    numMUs = size(spikeTrains,2);
    
    if normalization==1 % rec curve witwith generic curve  
                updatedTwitchFactor = zeros([1,size(spikeTrains,2)]);
%                 figure;
        r=[0,0.4,0.75,1];
       f=[0,0.2,0.5,1]; %normalized filtred force
        recCurve = fit(f',r','poly2');
% compare with
%Determine percentage of activation (duchateau and Enoka 2021) TA
%         f=0:0.01:1;
%             hold on,plot(f,recCurve(f));
% m = 0.0387;
% proportion = ((log(100*f)/m) - 2.1178)./100;
% 
% hold on, plot (f,proportion)
% 
% m = 0.0446;
% proportion = ((log(100*f)/m) - 2.1178)./100;
% hold on, plot (f,proportion)
% 
% rec_sol_MUs = [8,5,3,4,2,2,5,6,4,3];
% rec_sol_cumMUs = cumsum(rec_sol_MUs)/sum(rec_sol_MUs);
% rec_sol_activation = 0.025:.10:1;
% hold on, plot(rec_sol_cumMUs,rec_sol_activation)
% xlabel('Torque % MVC')
% ylabel('MU rec(0-1)')
% 
% legend('Old generic Curve', 'TA curve', 'Soleus', 'Soleus, (literature)') 
%     switch muscle
%         case 'TA',    totalMUs = 350;
%         case 'PERTERT',  totalMUs = 200; %check
%         case 'GASlat',  totalMUs = 500;
%         case 'GASmed',  totalMUs = 400;
%         case 'PERLONGUS',  totalMUs = 300; %check
%         case 'SOL',  totalMUs = 900;
%     end

    normFactor = 1./numMUs; % this goes together with recCurve
    elseif normalization ==2  % rec cruve according to Duchateau and Enoka 2020 Enokja https://journals.physiology.org/doi/epdf/10.1152/japplphysiol.00290.2021
        % and Oye et al.
        %numMUs=ones(size(numMUs)); % test get rid of it
        
         m = 0.045; % from duchateau and enoka 2022, so far I don't know models for other muscles , so going for the TA as generic 
        % WRONG : I used a scaled version of rafa for m //
%         switch muscle   obtained with script finalPlots, based on TA and SOL refereces below: 
%             case 'TA',    m = 0.0387;  totalMUs = 350; % Enoka to reach 1 at 50% recruited
%             case 'PERTERT', m = 0.0386; totalMUs = 200;
%             case 'GASlat',  m = 0.0366; totalMUs = 500;
%             case 'GASmed',  m = 0.0371; totalMUs = 400;
%             case 'PERLONGUS', m = 0.0396; totalMUs = 300;
%             case 'SOL',  m = 0.0446; totalMUs = 900; % oye et al to reach 1 of act at 95% recruited
%         end
        normFactor = 1./numMUs;
        %normFactor = (totalMUs-numMUs)./totalMUs;
        updatedTwitchFactor = zeros([1,size(spikeTrains,2)]);
        
        
    elseif normalization ==3  %only numMUs
        normFactor = 1./numMUs;
    elseif normalization ==4  %only numMUs
        normFactor = 1./numMUs;
        Ap = atanh(Ap);
        
        
        
    end





    %%
    for n = 3:length(spikeTrains)
        
        if sum(spikeTrains(n,:))>0
            n_spike=n_spike+1;
            t0 = (n-1)*T; % as n starts with 1
            firstSpike=1;
            t0s(n_spike,spikeTrains(n,:)>0)=t0;
        end
                
%         if n_spike >= 554  % just for debugging 
%             neuralMatrix(553,23)'
%         end
        if firstSpike % to ensure that only after the first spike compute the models 
            if strcmp(twitchModel,'Raikova') % round asures not taking very small negative values instead of 0 when the spike happens
                neuralMatrix(1:n_spike,:) = (p.*round(t(n)-t0s(1:n_spike,:),8)'.^m .*exp(-k.*round(t(n)-t0s(1:n_spike,:),8)'))';
                neuralActivation(n,:) = sum(neuralMatrix(1:n_spike,:),1,'omitnan');
            elseif strcmp(twitchModel,'Fuglevand')
                if normalization==0   % test og normalization
                    neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
                        + (((Ap*T)./tc).*exp(1-T./(tc)))'.*spikeTrains(n-1,:);
                elseif normalization ==1
                    % updateTwitchFactor updates the twitch only when a new
                    % a new twitch appears 
                    updatedTwitchFactor(logical(spikeTrains(n-1,:))) = 1.5*recCurve(abs(force(n)))+0.5;
                    neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
                        + ((((normFactor(i))*Ap*T)./tc).*exp(1-T./(tc)))'.*updatedTwitchFactor.*spikeTrains(n-1,:);
% practically the same the one above waits for next twitch 
                    %                     neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
%                         + ((((normFactor(i)*recCurve(force(n)))*Ap*T)./tc).*exp(1-T./(tc)))'.*excitation(n-1,:);
                 elseif normalization ==2
                    %normalizing according to recruitmentc coding
                    updatedTwitchFactor(logical(spikeTrains(n-1,:))) = ((log(100*abs(force(n-1)))/m) - 2.1178)./100;
                    if  updatedTwitchFactor(logical(spikeTrains(n-1,:)))>=1 %saturating to 1 
                        updatedTwitchFactor(logical(spikeTrains(n-1,:)))=1;
                    end
                    
                    %normalizing according to rate coding this might be
                    %wrong 
%                     if ~isnan(dischargeRate(n-1))
%                     updatedTwitchFactor(logical(excitation(n-1,:)))= dischargeRate(n-1).*updatedTwitchFactor(logical(excitation(n-1,:)));
%                     end
% if sum(~isnan(dischargeRate(n-1,find(logical(excitation(n-1,:))))))>=1
%     updatedTwitchFactor(logical(excitation(n-1,:)))= dischargeRate(n-1,find(logical(excitation(n-1,:)))).*...
%         updatedTwitchFactor(logical(excitation(n-1,:)));
% end
                    
                    neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
                        + ((((normFactor(i))*Ap*T)./tc).*exp(1-T./(tc)))'.*updatedTwitchFactor.*spikeTrains(n-1,:); 
                elseif normalization >=3
                    neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
                        + (((((normFactor(i)))*abs(force(n)-force(n-1))*Ap*T)./tc).*exp(1-T./(tc)))'.*spikeTrains(n-1,:);
                end
            elseif strcmp(twitchModel,'CEINMS')
                neuralActivation_(n,:) =  alpha_*spikeTrains(n,:) - beta1_*neuralActivation_(n-1,:)  - ...
                    beta2_*neuralActivation_(n-2,:);
                neuralActivation_(n,:) = activationScale_*(exp(shapeFactor_*neuralActivation_(n,:)) - 1) / (exp(shapeFactor_) - 1);
            elseif strcmp(twitchModel,'Piecewise') % does not work :(
                [rLessTc,cLessTc] = find(t(n)-t0s(1:n_spike,:) < tc);
                [rInTcThr,cInTcThr] = find(((t(n)-t0s(1:n_spike,:) > tc) .* (t(n)-t0s(1:n_spike,:)< 2*(thr-tc)+tc))>0);
                if n_spike>=2
                    [rMoreThr,cMoreThr] = find(t(n)-t0s(1:n_spike,:)> 2*(thr-tc)+tc);
                    neuralMatrix(rMoreThr,cMoreThr) =  0;
                end
                neuralMatrix(rLessTc,cLessTc) = 0.5*(Ap-Ap*cos(pi./(tc)*(t(n)-t0s(rLessTc,cLessTc))));
                neuralMatrix(rInTcThr,cInTcThr) = ...
                    0.5*(Ap-Ap*cos(pi./(2*(thr-tc))*(t(n)-t0s(rInTcThr,cInTcThr)+2*(thr-3*tc/2))));
                
                neuralActivation(n,:) = sum(neuralMatrix(1:n_spike,:),1,'omitnan');
            end
            if tetanusModel
                neuralActivationTransf(n,:) = normSat*(1-1*exp(-tetParam*(neuralActivation(n,:))))./(1+1*exp(-tetParam*(neuralActivation(n,:))));
            end
        end
    end
%     size(excitation,2)
%                     totalMUs*recCurve(max(force))
%                     totalMUs*recCurve(max(force))/size(excitation,2)
    % hold on, plot(t,neuralActivation)
    if tetanusModel
        neuralActivation = neuralActivationTransf;
    end
end
