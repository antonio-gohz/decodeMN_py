function plotMUAP(MUAP,varargin)
% Plot 2D motor unit action potentials

channelConfig = [{1:8},{9:16},{17:24},{25:32},{33:40},{41:48},{49:56},{57:64}];
fs=2048;
scaleX = size(MUAP,2);
scaleY = 2*max(abs(MUAP(:)));
%color = [0.7,0.7,0.7];
color = {'#4477AA'};
% Process optional input parameters using a loop and switch case
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i+1};
    switch param
        case 'fs'
            fs = value;
        case 'scaleX'
            scaleX = value;
        case 'scaleY'
            scaleY = value;
        case 'channelConfig'
            channelConfig = value;
        case 'Color'
            color = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end
%MUAP=MUAP/scaleY;
%TF = isoutlier(max(abs(MUAP),[],2));
%TF = max(abs(MUAP),[],2)<L;
offset=0;
p=offset;


for i=1:length(channelConfig)
    for CH = channelConfig{i}
%         if ~TF(CH)
%             color='b';
%         else
%             color='r';
%         end
        hold on,plot([0:size(MUAP,2)-1]/fs+i*scaleX/fs,MUAP(CH,:)-p*scaleY,'Color',color);
        p=p+1;
    end
    p=0;
end
cols =1:size(channelConfig,2);
rows = 1:max(cellfun(@(x)size(x,2), channelConfig));
xticks(cols.*scaleX/fs+scaleX/(2*fs))
yticks(flip((rows-1).*scaleY*-1))
xticklabels(cols);
yticklabels(rows);
xlabel('Columns')
ylabel('Rows')