function hAX = multipleAxes(hAX,axesLimits,xLabels,yLabel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% figure;
% hAX=axes;                 % first axes, save handle
pos=get(hAX,'position');   % get the position vector
pos1=pos(2);              % save the original bottom position
pos(2)=pos(2)+2*pos1; pos(4)=pos(4)-2*pos1;  % raise bottom/reduce height->same overall upper position
set(hAX,'position',pos)   % and resize first axes (raise it)
pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
pos1=pos(2);              % save the original bottom position
pos(2)=pos(2)+pos1; pos(4)=0;  % raise bottom/reduce height->same overall upper position
hAX(2)=axes('position',pos,'color','none');  % and create the second
pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
hAX(3)=axes('position',pos,'color','none');  % and create the second


xlabel(hAX(1),xLabels{1})
ylabel(hAX(1),yLabel)
linkaxes(hAX)

xlabel(hAX(2),xLabels{2})
xticklabels(hAX(2),round(linspace(axesLimits(1,1), axesLimits(1,2),length(hAX(1).XTick) ),1))
set(hAX(2),'FontSize',10,'linewidth',1,'box','off')

xlabel(hAX(3), xLabels{3})
xticklabels(hAX(3),round(linspace(axesLimits(2,1), axesLimits(2,2),length(hAX(1).XTick) )))
set(hAX(3),'FontSize',10,'linewidth',1,'box','off')

% set(haX(2),'xcolor','r','ycolor','r')
end

