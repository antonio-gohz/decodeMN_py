function [h,ax] = updatePlot(data,h,ax,xlimits,formats)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
 for i=1:size(h,2)
     if ~isempty(data{i,2})
     set(h(:,i), 'XData',data{i,1},...
         {'YData'},num2cell(data{i,2},2)')
     drawnow limitrate;
     end
 end
 
if nargin>3
    set(ax,{'Xlim'},{[xlimits(1) xlimits(2)]})
end
end
