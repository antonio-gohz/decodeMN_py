function Fig = scatter_distribution_figure(datax, datay, cats, varargin)

% Function to create a scatter plot showing the individual data points
% and the distributions on the X and Y axis. 
%
% INPUT
% scatter_distribution_figure(datax, datay, cats)
% datax:    N x 1 vector containing x-values of the data to be plotted
% datay:    N x 1 vector containing y-values of the data to be plotted 
% cats:     N x 1 cell array of character vectors representing the 
%           corresponding groups of the data points
%
% OPTIONAL INPUT
% scatter_plot(..., 'PARAM1', val1, 'PARAM2', val2, ...)
% 
%     'Colors'       Color of the jitter plots, num of cells should 
%                    correspond to categories{[r,g,b], [r,g,b], ...}. 
%                    Defaults to copper colormap 
%     'Markersize'   Markersize of data points.
%                    Defaults to 100. 
%     'Trendline'    on or off. 
%                    Defaults to off.  
%     'XLim'         Xlim [min max]
%                    Defaults to standard lims. 
%     'YLim'         Ylim [min max]
%                    Defaults to standard lims. 
%     'XLabel'       Xlabel {str}. 
%                    Defaults to empty string. 
%     'YLabel'       Ylabel {str}. 
%                    Defaults to empty string.
%     'DistType'     Kernel or Gaussian. 
%                    Defaults to Kernel.
%
% Copyright (c) 2022, Eline Zwijgers, Sint Maartenskliniek, 
%                     e.zwijgers@maartenskliniek.nl

% Input erros 
if nargin < 3
   error('Not enough input arguments'); 
end

if length(datax) ~= length(datay)
   error('Datax and datay vector should be the same length'); 
end

if length(datax) ~= length(datay)
   error('Data and category vector should be the same length'); 
end

% Category settings
cats = categorical(cats);
catnames = (unique(cats,'stable')); 

% Default plot settings 
cols = mat2cell(colormap(lines(length(catnames))), ones(1,length(catnames)), 3)';
marker_size = 100; 
trendline = 'off'; 
plot_type = 'Kernel'; 

% Optional plot settings 
for j = 1:length(varargin)
   if isequal(varargin{j},'Colors')
      cols = varargin{j+1}; 
      if length(cols) ~= length(catnames) 
          warning('Number of colors is not equal to the number of categories. The default colors are used') 
          cols = mat2cell(colormap(lines(length(catnames))), ones(1,length(catnames)), 3)';
      end  
   end
   if isequal(varargin{j},'Markersize')
      marker_size = varargin{j+1}; 
   end
   if isequal(varargin{j},'Trendline')
      trendline = varargin{j+1}; 
   end 
   if isequal(varargin{j},'YLim')
      y_lim = varargin{j+1}; 
   end 
   if isequal(varargin{j},'XLim')
      x_lim = varargin{j+1}; 
   end 
   if isequal(varargin{j},'XLabel')
      x_label = varargin{j+1}; 
   end 
   if isequal(varargin{j},'YLabel')
      y_label = varargin{j+1}; 
   end 
   if isequal(varargin{j},'PlotType') 
        if ~strcmp(dist_type,'Kernel') && ~strcmp(dist_type,'Gaussian')
            warning('Unknown DistType. Choose between "Kernel" and "Gaussian". The default distribution (Kernal) is used.')
            dist_type = 'Kernel';
        end
   end 
end

% Scatter plot 
Fig.scat = subplot(4,5,[6 7 8 9 11 12 13 14 16 17 18 19]); 
hold on 
for n = 1:length(catnames)
    thisCat = catnames(n);
    thisDatax = datax(cats == thisCat);
    thisDatay = datay(cats == thisCat);
    Scatter(thisDatax,thisDatay, cell2mat(cols(n)), marker_size, trendline);
end

if exist('x_lim','var') 
    set(Fig.scat,'XLim',x_lim) 
end 

if exist('y_lim','var') 
    set(Fig.scat,'YLim',y_lim) 
end 

if exist('x_label','var')
    set(Fig.scat.XAxis.Label,'String',x_label)
end 

if exist('y_label','var')
    set(Fig.scat.YAxis.Label,'String',y_label)
end 

switch trendline
    case 'off'
        legend(catnames,'Location','northwest')
    case 'on'
        legend(Fig.scat.Children(length(catnames)*2:-2:1,1), ...
        catnames,'Location','northwest')   
end 

% Distribution plot y-axis 
Fig.disty = subplot(4,5,[10 15 20]); 
axis tight
hold on 
for n = 1:length(catnames)
    thisCat = catnames(n);
    thisData = datay(cats == thisCat);
    Distribution(thisData, cell2mat(cols(n)),'y', plot_type); 
end

% Format axis of distribution plot y-axis 
if exist('y_lim','var')
    set(Fig.disty,'YLim',y_lim)
else
    yLimits = get(Fig.disty,'YLim');
    set(Fig.scat,'ylim',yLimits)
end

set(Fig.disty,'xtick',[],'ytick',[],'Xcolor','none','Ycolor','none')

% Distribution plot x-axis 
Fig.distx = subplot(4,5,[1 2 3 4]); 
axis tight
hold on 
for n = 1:length(catnames)
    thisCat = catnames(n);
    thisData = datax(cats == thisCat);
    Distribution(thisData, cell2mat(cols(n)),'x', plot_type); 
end

% Format axis of distribution plot x-axis
if exist('x_lim','var')
    set(Fig.distx,'XLim',x_lim)
else
    xLimits = get(Fig.distx,'XLim');
    set(Fig.scat,'xlim',xLimits)
end

set(Fig.distx,'xtick',[],'ytick',[],'Xcolor','none','Ycolor','none')


% Functions 
function Scatter(datax, datay, col, marker_size, trendline)
    scatter(datax, datay, marker_size, 'filled','MarkerEdgeColor','None','MarkerFaceColor',col,'MarkerFaceAlpha',1);

    if contains(trendline,'on')
        lf = fitlm(datax,datay);
        fit = @(x) lf.Coefficients{2,1}*x + lf.Coefficients{1,1}; 
        lim = [min(datax)-0.05 max(datax)+0.05];
        p = plot(lim,fit(lim),'color',col,'LineWidth',2);
        p.Color(4) = 0.6;
    end 
end 

function Distribution(data, col, dir, plot_type)
    mean_data = mean(data, 'omitnan'); 
    std_data = std(data, 'omitnan'); 
    xnormdis = (-3*std_data+mean_data:0.001:3*std_data+mean_data); 
    y_norm = normpdf(xnormdis,mean_data,std_data);
    
    switch plot_type 
        case 'Gaussian'
            switch dir 
                case 'y'
                    plot(y_norm/10,xnormdis,'Color',col,'LineWidth',1.5)
                    y_norm(1) = 0; 
                    y_norm(end) = 0; 
                    patch(y_norm/10,xnormdis,col,'FaceAlpha',.4,'EdgeColor','none')
                case 'x'
                    plot(xnormdis,y_norm/10,'Color',col,'LineWidth',1.5)
                    y_norm(1) = 0; 
                    y_norm(end) = 0; 
                    patch(xnormdis,y_norm/10,col,'FaceAlpha',.4,'EdgeColor','none')
            end 
 
        case 'Kernel'
            [f, xi] = ksdensity(data); 
            switch dir
                case 'y'
                    plot(f,xi,'Color',col,'LineWidth',1.5)
                    patch(f,xi,col,'FaceAlpha',0.4,'EdgeColor','none')
                case 'x'
                    plot(xi,f,'Color',col,'LineWidth',1.5)
                    patch(xi,f,col,'FaceAlpha',0.4,'EdgeColor','none')
            end 
    end          
            
end 
    
end 