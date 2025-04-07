function name = check_overwriting(name)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
counter=1;
% save data
while 1
    if exist([name,'.mat'],'file')==2
        name = [name,'_',num2str(counter)];
        counter = counter+1;
    else
        break;
    end
end
end