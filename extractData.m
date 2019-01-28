function [  ] = extractData(output, analysisplan)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
cd(output)
[num,~,raw] = xlsread(analysisplan,'experiments');
listnums = find(num(:,5)==0)+1;

for filenum = 1:length(listnums)
    curfile = num2str(raw{listnums(filenum),1});
    load([curfile '_rawdata.mat']);
    load([curfile '_metadata.mat']);
data = zeros(size(rawdata,1),1);
for i = 1:size(rawdata,1)
    % series number is in last column
    data(i) = metadata.time(rawdata(i,end),rawdata(i,3));
end
% data = [metadata.um2px*rawdata(:,1:2), rawdata(:,[3,4,end]), data]; % without [collagen] or treatment info
data = [metadata.um2px*rawdata(:,1:2), rawdata(:,[3,4,end-2:end]), data]; % with [collagen] and treatment info 
save([curfile '_data'],'data');
xlswrite(analysisplan,1,'experiments',['H' num2str(listnums(filenum))]);
end
if(isempty(listnums))
    fprintf('Data for experiments in the list were already extracted\n')
end
end