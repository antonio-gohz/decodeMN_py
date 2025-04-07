function [hg,pd,data_pdfz] = histPDF(data,binEdges,color)
%plotHistPDF plots histogram and PDFs
%   Detailed explanation goes here

hg = histogram(data,binEdges,'LineStyle','none','Normalization','pdf','FaceColor', color);
%GMModel = fitgmdist(data,2);
pd = fitdist(data,'Kernel');%,'Kernel','epanechnikov'); %,'Kernel','epanechnikov');
x_data = linspace(hg.BinLimits(1),hg.BinLimits(2),50); % x axis for histogram
x_data_pdf = linspace(min(x_data),max(x_data),1000);
data_pdfz = pdf(pd,x_data_pdf');
%data_pdfz = data_pdfz*max(hg.Values)/max(data_pdfz); % for visualisation matchign maximum values

plot(x_data_pdf,data_pdfz,'LineWidth', 1.5, 'HandleVisibility', 'off', 'Color', color);
data_pdfz=[x_data_pdf',data_pdfz];
end

