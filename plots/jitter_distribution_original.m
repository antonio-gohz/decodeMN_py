function Fig = jitter_distribution_pairs(data, cats, varargin)

% Function to create a figure showing the individual data points, mean, 
% standard deviation, and distribution of multiple categories. 
%
% INPUT
% jitter_distribution_figure(data, cats)
% data:     N x 1 vector containing the data points to be plotted
% cats:     N x 1 cell array of character vectors representing the 
%           corresponding groups of the data points
%
% OPTIONAL INPUT
% jitter_plot(..., 'PARAM1', val1, 'PARAM2', val2, ...)
% 
%     'Colors'       Color of the jitter plots, num of cells should 
%                    correspond to categories{[r,g,b], [r,g,b], ...}. 
%                    Defaults to lines colormap. 
%     'Markersize'   Markersize of data points.
%                    Defaults to 200.   
%     'Linewidth'    Linewidth of the mean and std plot.
%                    Defaults to 3.  
%     'Capsize'      Capsize of the mean and std plot. 
%                    Defaults to 15.
%     'YLim'         Ylim [min max]
%                    Defaults to standard lims. 
%     'YLabel'       Ylabel {str}. 
%                    Defaults to empty string.
%     'DistType'     'Kernel' or 'Gaussian'. 
%                    Defaults to Kernel. 
%     'PlotType'     'Interal' or 'External'.
%                    The option 'Interal' plots the distribution to the
%                    right of the mean and errorbar. 
%                    The option 'External' plots the distribution to the
%                    right of the figure outline.   
%                    Defaults to External. 
%
% Copyright (c) 2022, Eline Zwijgers, Sint Maartenskliniek, 
%                     e.zwijgers@maartenskliniek.nl

% Input erros 
if nargin < 2
   error('Not enough input arguments'); 
end

if length(data) ~= length(cats)
   error('Data and category vector should be the same length'); 
end

% Category settings
cats = categorical(cats);
catnames = (unique(cats,'stable'));    
catnames_labels = char(catnames);

% Default plot settings 
cols = mat2cell(colormap(lines(length(catnames))), ones(1,length(catnames)), 3)';
marker_size = 200; 
line_width = 3;
cap_size = 15;
dist_type = 'Kernel'; 
plot_type = 'External';

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
   if isequal(varargin{j},'Linewidth')
      line_width = varargin{j+1}; 
   end
   if isequal(varargin{j},'Capsize')
      cap_size = varargin{j+1}; 
   end
   if isequal(varargin{j},'YLim')
      y_lim = varargin{j+1}; 
   end 
   if isequal(varargin{j},'YLabel')
      y_label = varargin{j+1}; 
   end 
   if isequal(varargin{j},'DistType') 
      dist_type = varargin{j+1};
        if ~strcmp(dist_type,'Kernel') && ~strcmp(dist_type,'Gaussian')
            warning('Unknown DistType. Choose between "Kernel" and "Gaussian". The default distribution (Kernal) is used.')
            dist_type = 'Kernel'; 
        end 
   end 
   if isequal(varargin{j},'PlotType') 
      plot_type = varargin{j+1};
        if ~strcmp(plot_type,'Internal') && ~strcmp(plot_type,'External')
            warning('Unknown PlotType. Choose between "External" and "Internal". The default plot setting (External) is used.')
            plot_type = 'External'; 
        end 
   end 
end

% Plotting 
switch plot_type 
    case 'External'
        % Jitter plot 
        Fig.jit = subplot(6,4,[1 2 3 5 6 7 9 10 11 13 14 15 17 18 19 21 22 23]); 
        hold on 

        for n = 1:length(catnames)
            thisCat = catnames(n);
            thisData = data(cats == thisCat);
            Jitter(thisData, n, cell2mat(cols(n)), marker_size, line_width, cap_size);
        end
            
        % Format axis of jitter plot 
        set(Fig.jit, 'XTick', (1:length(catnames))+0.1, 'XTickLabels', catnames_labels);
        set(Fig.jit,'XLim',[0.5 length(catnames)+0.6])  

        if exist('y_lim','var') 
            set(Fig.jit,'YLim',y_lim) 
        end 

        if exist('y_label','var')
            set(Fig.jit.YAxis.Label,'String',y_label)
        end 

        % Distribution plot 
        Fig.dist = subplot(6,4,[4 8 12 16 20 24]); 
        axis tight
        hold on 

        for n = 1:length(catnames)
            thisCat = catnames(n);
            thisData = data(cats == thisCat);
            Distribution(thisData, cell2mat(cols(n)), line_width, dist_type); 
        end
        
        % Format axis of distribution plot 
        if exist('y_lim','var')
            set(Fig.dist,'YLim',y_lim)
        else
            yLimits = get(Fig.dist,'YLim');
            set(Fig.jit,'ylim',yLimits)
        end
        set(Fig.dist,'xtick',[],'ytick',[],'Xcolor','none','Ycolor','none')

    case 'Internal'
        % Jitter plot combined with distribution plot  
        Fig.jit = subplot(1,1,1); 
        hold on 

        for n = 1:length(catnames)
            thisCat = catnames(n);
            thisData = data(cats == thisCat);
            
            % Scale of distribution 
            [ydens, ~] = ksdensity(data);
            switch dist_type
                case 'Kernel'
                    scale = 0.20/max(ydens); 
                case 'Gaussian'
                    scale = 0.1/max(ydens); 
            end 

            Jitter_distribution(thisData, n, cell2mat(cols(n)), marker_size, line_width, cap_size, dist_type, scale);

        end
        
        % Format axis 
        set(Fig.jit, 'XTick', (1:length(catnames))+0.1, 'XTickLabels', catnames_labels);

        set(Fig.jit,'XLim',[0.5 length(catnames)+0.8]) 

        if exist('y_lim','var') 
            set(Fig.jit,'YLim',y_lim) 
        end 

        if exist('y_label','var')
            set(Fig.jit.YAxis.Label,'String',y_label)
        end 
              
end 
        
% Functions 

function Jitter(data, pos, col, marker_size, line_width, cap_size)
    
    [density, value] = ksdensity(data);          
    density = density(value >= min(data) & value <= max(data));
    value = value(value >= min(data) & value <= max(data));
    value(1) = min(data);
    value(end) = max(data);

    width = 0.05/max(density);
    jitterstrength = interp1(value, density*width, data);
    jit = 2*(rand(size(data))-0.5);

    scatter(pos + jit.*jitterstrength, data, marker_size, 'filled','MarkerEdgeColor','None','MarkerFaceColor',col,'MarkerFaceAlpha',0.6);    
    errorbar(pos+0.20,mean(data),std(data),'Color','k','LineWidth',line_width,'CapSize',cap_size);
    scatter(pos+0.20, mean(data),marker_size,'filled','MarkerFaceColor','k','MarkerEdgeColor', 'None');

end 

function Distribution(data, col, line_width, plot_type)
    
     switch plot_type 
        case 'Gaussian'
            mean_data = mean(data, 'omitnan'); 
            std_data = std(data, 'omitnan'); 
            xnormdis = (-3*std_data+mean_data:0.001:3*std_data+mean_data); 
            y_norm = normpdf(xnormdis,mean_data,std_data);
            plot(y_norm,xnormdis,'Color',col,'LineWidth',0.7*line_width)
            y_norm(1) = 0; 
            y_norm(end) = 0; 
            patch(y_norm,xnormdis,col,'FaceAlpha',0.4,'EdgeColor','none')
         case 'Kernel'
            [f, xi] = ksdensity(data); 
            plot(f,xi,'Color',col,'LineWidth',0.7*line_width)
            patch(f,xi,col,'FaceAlpha',0.4,'EdgeColor','none')
     end
     
end 

end 


function Jitter_distribution(data, pos, col, marker_size, line_width, cap_size, plot_type, scale)
    [density, value] = ksdensity(data);          
    density = density(value >= min(data) & value <= max(data));
    value = value(value >= min(data) & value <= max(data));
    value(1) = min(data);
    value(end) = max(data);

    width = 0.05/max(density);
    jitterstrength = interp1(value, density*width, data);
    jit = 2*(rand(size(data))-0.5);
    
    scatter(pos + jit.*jitterstrength, data, marker_size, 'filled','MarkerEdgeColor','None','MarkerFaceColor',col,'MarkerFaceAlpha',0.6);

    switch plot_type 
        case 'Gaussian'
            mean_data = mean(data, 'omitnan'); 
            std_data = std(data, 'omitnan'); 
            xnormdis = (-3*std_data+mean_data:0.001:3*std_data+mean_data); 
            y_norm = normpdf(xnormdis,mean_data,std_data);
            y_norm(1) = 0; 
            y_norm(end) = 0; 
            y_norm = y_norm*scale; 
            plot(y_norm+pos+0.20,xnormdis,'Color',col,'LineWidth',0.7*line_width)
            patch(y_norm+pos+0.20,xnormdis,col,'FaceAlpha',0.4,'EdgeColor','none')
        case 'Kernel'
            [f, xi] = ksdensity(data); 
            f = f*scale; 
            min_f = min(f); 
            offset = pos+0.2-min_f; 
            plot(f+offset,xi,'Color',col,'LineWidth',0.7*line_width)
            patch(f+offset,xi,col,'FaceAlpha',0.4,'EdgeColor','none')
    end 
   
    errorbar(pos+0.20,mean(data),std(data),'Color','k','LineWidth',line_width,'CapSize',cap_size);
    scatter(pos+0.20, mean(data),marker_size,'filled','MarkerFaceColor','k','MarkerEdgeColor', 'None');

end 



