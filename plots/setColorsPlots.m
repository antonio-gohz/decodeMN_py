function [colorsRGB,colorsHEX] = setColorsPlots(colorset)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

% default 
if nargin<1
colorset = 'BrightQualitative';
end

switch colorset
    case 'BrightQualitative'
    colorsHEX = {'#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'}; % Bright qualitative
end
colorsRGB =  hex2rgb(colorsHEX);