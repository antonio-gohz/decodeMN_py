function pairedCombs = subgroupPairs_v2(N,W)
% subgroupPairs_v1.m
%test better combos
% max number of unique combos with non repetable elements is
% factorial(N)/(factorial(W)*factorial(N-W))/2  in our case W is N/2
% factorial(N)/(factorial(N/2)*factorial(N/2))/2 
% for N=8, maxCombs=35;

%N=6;
% todo add W and if its not declare 
if nargin < 2 || isempty(W)
    W = floor(N/2);
end
pairedCombs{1}(1,:) = 1:W;
pairedCombs{1}(2,:) = W+1: 2*W;
maxCombs=2000;
% for generalization (although in this research we don't take data with
% less than 8 MUs i.e., maxCombs = 35; 
if maxCombs> factorial(N)/(factorial(W)*factorial(W))/2 
    maxCombs = factorial(N)/(factorial(W)*factorial(W))/2 ;
end
k=2;
flag_repetead=0;
rep=1;
while k<=maxCombs && rep<=maxCombs
    combs = randperm(N);
    for i=1:length(pairedCombs)
        % if current combs already exists in the pairedCombs
        if isequal(sort(combs(1:W)),sort(pairedCombs{i}(1,:))) || ... % combs 1 vs pairedCombs 1 
                isequal(sort(combs(1:W)),sort(pairedCombs{i}(2,:))) || ... % combs 1 vs pairedCombs 2 
                isequal(sort(combs(W+1:2*W)),sort(pairedCombs{i}(1,:))) || ... % combs 2 vs pairedCombs 1 
                isequal(sort(combs(W+1:2*W)),sort(pairedCombs{i}(2,:)))% combs 2 vs pairedCombs 2
            flag_repetead=1;
            break;
        end
    end
    if ~flag_repetead       
        pairedCombs{k}(1,:)=sort(combs(1:W));
        pairedCombs{k}(2,:)=sort(combs(W+1:2*W));
        k=k+1;
    else
        flag_repetead=0; %reset flag for next iteration
        rep=rep+1;
    end
end
