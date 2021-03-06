function [  ] = extractFeatures(output, analysisplan)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Update: with the 2 extra columns ([collagen] and treatment), I changed
%         the time column index from 6 to 8 (lines 29, 33, 39) (P.A. 6/16/18)
%
% Update: found bug on lines 21 and 23 (column 5 had changed to column 7)
% (6/19/2018)

cd(output)

[num,~,raw] = xlsread(analysisplan,'experiments');
filenums = find(num(:,6)==0)+1;

% figure(100)

% downsampling_factor = input('Please input downsampling factor: ');

for curfilenum = 1:length(filenums)
    curfile = [num2str(raw{filenums(curfilenum),1}) '_data'];
    load(curfile);
    SPEED = [];
    DELTA_THETA = [];
    AV_SPEED = [];
    DIST = [];
    stages = unique(data(:,7));
    for i=1:length(stages)
        stageind = find(data(:,7)==stages(i));
        stagedata = data(stageind,:);
        cells = unique(stagedata(:,4));
        for j=1:length(cells)
            cellind = find(stagedata(:,4)==cells(j));
            celldata = stagedata(cellind,:);
%             celldata = downsample(celldata, downsampling_factor);
%             plot(celldata(:,1),celldata(:,2)); grid on; axis square;
%             hold on
            diffcelldata = diff(celldata);
            speedvec = sqrt(diffcelldata(:,1).^2+diffcelldata(:,2).^2)./diffcelldata(:,8);
            speedvec = [nan; speedvec];
            SPEED = [SPEED; speedvec];
            distvec = sqrt(celldata(:,1).^2+celldata(:,2).^2);
            DIST = [DIST; distvec];
            theta = atan2d(diffcelldata(:,1),diffcelldata(:,2)); %absolute direction of the cell
            persistence = zeros(size(celldata,1)-2,1);
            for k = 1:(size(celldata,1)-2)
                ind1 = k;
                ind3 = k+2;
%                 persistence(k) = sqrt((celldata(ind3,1)-celldata(ind1,1))^2+(celldata(ind3,2)-celldata(ind1,2))^2)/(celldata(ind3,8)-celldata(ind1,8)); % magnitude of average velocity (total displacement) between point 1 and 3, 2 and 4 etc..
                ind2 = k+1;
                dist12 = sqrt((celldata(ind2,1)-celldata(ind1,1))^2+(celldata(ind2,2)-celldata(ind1,2))^2); % distance between points 1 and 2
                dist23 = sqrt((celldata(ind3,1)-celldata(ind2,1))^2+(celldata(ind3,2)-celldata(ind2,2))^2); % distance between points 2 and 3
                persistence(k) = (dist12 + dist23)/(celldata(ind3,8)-celldata(ind1,8)); % average speed (total distance traveled / time) between point 1 and 3, 2 and 4 etc..
            end
            persistence = [nan; nan; persistence];
%             theta = [nan; theta];
            delta_theta = [nan; nan; diff(theta)];
            % force values to be between -180 and 180
            delta_theta(delta_theta>180) = delta_theta(delta_theta>180)-360;
            delta_theta(delta_theta<-180) = delta_theta(delta_theta<-180)+360;
            delta_theta = abs(delta_theta);
            AV_SPEED = [AV_SPEED; persistence];
            DELTA_THETA = [DELTA_THETA; delta_theta];
       end
    end
features = [data(:,[1,2,8]), SPEED, DELTA_THETA, AV_SPEED, data(:,5:7), DIST];
save([curfile(1:end-4) 'features'],'features');
xlswrite(analysisplan,1,'experiments',['I' num2str(filenums(curfilenum))]);

if(isempty(listnums))
    fprintf('Features for experiments in the list were already extracted\n')
end

end

% hold off