function hAY = multipleYAxes(hAY,axesLimits,yLabels )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% figure;
% hAX=axes;                 % first axes, save handle
pos=get(hAY,'position');   % get the position vector
pos1=pos(1);              % save the original bottom position
pos(1)=pos(1)+0.05; pos(3)=pos(3)-0.05;  % raise bottom/reduce height->same overall upper position
set(hAY,'position',pos)   % and resize first axes (raise it)
pos(1)=pos1; pos(3)=0; % reset bottom to original and small height
hAY(2)=axes('position',pos,'color','none');  % and create the second
pos(1)=pos1-0.05; pos(3)=0; % reset bottom to original and small height
hAY(3)=axes('position',pos,'color','none');  % and create the second

ylabel(hAY(1),yLabels{1})
linkaxes(hAY)

ylabel(hAY(2),yLabels{2})
yticklabels(hAY(2),round(linspace(axesLimits(1,1), axesLimits(1,2),length(hAY(1).YTick) ),1))
set(hAY(2),'FontSize',14,'linewidth',1,'box','off')

ylabel(hAY(3), yLabels{3})
yticklabels(hAY(3),round(linspace(axesLimits(2,1), axesLimits(2,2),length(hAY(1).YTick) )))
set(hAY(3),'FontSize',14,'linewidth',1,'box','off')


% set(haX(2),'xcolor','r','ycolor','r')
end