function [  ] = extractFeatures( output )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Update: with the 2 extra columns ([collagen] and treatment), I changed
%         the time column index from 6 to 8 (lines 29, 33, 39) (P.A. 6/16/18)
%
% Update: found bug on ines 21 and 23 (column 5 had changed to column 7)
% (6/19/2018)

% output = '\\PHYS34212\MigrationData\MigrationData\Migration1\Output';
cd(output)
% analysisplan = '\\PHYS34212\MigrationData\MigrationData\Migration1\Code files v 2\analysisplan.xlsx';
analysisplan = '\\PHYS34212\MigrationData\MigrationData\Migration1\Code files v 2\analysisplan.xls';

[num,~,raw] = xlsread(analysisplan,'experiments');
filenums = find(num(:,6)==0)+1;

for curfilenum = 1:length(filenums)
    curfile = [num2str(raw{filenums(curfilenum),1}) '_data'];
    load(curfile);
    SPEED = [];
    THETA = [];
    PERSISTENCE = [];
    stages = unique(data(:,7));
    for i=1:length(stages)
        stageind = find(data(:,7)==stages(i));
        stagedata = data(stageind,:);
        cells = unique(stagedata(:,4));
        for j=1:length(cells)
            cellind = find(stagedata(:,4)==cells(j));
            celldata = stagedata(cellind,:);
            plot(celldata(:,1),celldata(:,2)); grid on; axis square;
            diffcelldata = diff(celldata);
            speedvec = sqrt(diffcelldata(:,1).^2+diffcelldata(:,2).^2)./diffcelldata(:,8);
            speedvec = [nan; speedvec];
            SPEED = [SPEED; speedvec];
            theta = atan2d(celldata(:,1),celldata(:,2));
            persistence = diff(theta)./diffcelldata(:,8);
            persistence = [nan; persistence];
            PERSISTENCE = [PERSISTENCE; persistence];
            THETA = [THETA; theta];
       end
    end
features = [data(:,[1,2,8]), SPEED, THETA, PERSISTENCE, data(:,5:7)];
save([curfile(1:end-4) 'features'],'features');
xlswrite(analysisplan,1,'experiments',['I' num2str(filenums(curfilenum))]);
end