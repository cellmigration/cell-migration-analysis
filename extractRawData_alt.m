function [  ] = extractRawData_alt ( trackpath , output)
%Reads tracking data from mat file(s), extracts x,y location of each
%tracked cell in each stage (series). The output is a matrix where each row
%is defined as following: 1)relative y-coordinate 2)relative x-coordinate
%3)frame number 4)cell ID 5)series number

%   trackpath is the path where the mat files are found and outputpath is
%   where the raw variable will be saved. Each tracking data file should
%   be named based on this template:'yyyymmddCOL123C123S12' (it is the nd2
%   file name + 'S01' that represents the stage (series) number.

% Update: now extracts raw data from tracking files for each row in the 
%         analysis plan. Series collagen concentration, treatment and
%         series number are in the last 3 columns of raw data matrix.
%
% Update: adapt for old PC3 data by uncommenting lines 43-44 and 77-78

% list of experiments
% analysisplan = '\\PHYS34212\MigrationData\MigrationData\Migration1\code files v 2 check\analysisplan';
analysisplan = '\\phys34212\migrationdata\MigrationData\Migration1\Repositories\cell-migration-analysis\analysisplan.xls';

[num,~,raw] = xlsread(analysisplan,'experiments');
listnums = find(num(:,3)==0)+1;

if(isempty(listnums))
    fprintf('Tracking data for experiments in the list were already extracted\n')
else

rawdata = [];
prevfolder = num2str(raw{listnums(1),1});


for numofseries = 1:length(listnums)
    % Change directory to the one containing series
    curfolder = num2str(raw{listnums(numofseries),1});
    if strcmp(curfolder,prevfolder)
        cd([trackpath '\' curfolder]);
        % Read series conditions (sample naming convention: S001)
        seriesname = raw{listnums(numofseries),3};
        seriesnumber = str2double(seriesname(2:end));
        
        % Uncomment for old PC3 data (sample naming convention: 20171107HAC1S47)
%         seriesname = raw{listnums(numofseries),3};
%         seriesnumber = str2double(seriesname(end-1:end));
        
        
        seriescollconc = raw{listnums(numofseries),4}; % collagen concentration
        seriestreat = raw{listnums(numofseries),5}; % experimental treatment (1 = control, 2 = HA)
        seriesfilename = [seriesname, '.mat'];
        % load series
        load(seriesfilename);
        cellnums = unique(pos(:,4));
        % Add series collagen concentration, treatment and number
        pos = [pos, seriescollconc*ones(size(pos,1),1),seriestreat*ones(size(pos,1),1), seriesnumber*ones(size(pos,1),1)];
        for curcellnum = 1:length(cellnums)
            curcellind = find(pos(:,4)==cellnums(curcellnum));
            curcellrawdata = pos(curcellind,:);
            % changing coordinates relative to origin of cell
            curcellrawdata(:,1) = -curcellrawdata(:,1)+curcellrawdata(1,1);
            curcellrawdata(:,2) = curcellrawdata(:,2)-curcellrawdata(1,2);
            % concatenating cell data to rawdata
            rawdata = [rawdata; curcellrawdata];
        end
        prevfolder = curfolder;
        
    else
        cd(output)
        save([prevfolder '_rawdata'],'rawdata');
        rawdata = [];
        cd([trackpath '\' curfolder]);
        
        % Read series conditions (sample naming convention: S001)
        seriesname = raw{listnums(numofseries),3};
        seriesnumber = str2double(seriesname(2:end));
        
        % Uncomment for old PC3 data (sample naming convention: 20171107HAC1S47)
%         seriesname = raw{listnums(numofseries),3};
%         seriesnumber = str2double(seriesname(end-1:end));
      
        seriescollconc = raw{listnums(numofseries),4}; % collagen concentration
        seriestreat = raw{listnums(numofseries),5}; % experimental treatment (1 = control, 2 = HA)
        seriesfilename = [seriesname, '.mat'];
        % load series
        load(seriesfilename);
        cellnums = unique(pos(:,4));
        % Add series collagen concentration, treatment and number
        pos = [pos, seriescollconc*ones(size(pos,1),1),seriestreat*ones(size(pos,1),1), seriesnumber*ones(size(pos,1),1)];
        for curcellnum = 1:length(cellnums)
            curcellind = find(pos(:,4)==cellnums(curcellnum));
            curcellrawdata = pos(curcellind,:);
            % changing coordinates relative to origin of cell
            curcellrawdata(:,1) = -curcellrawdata(:,1)+curcellrawdata(1,1);
            curcellrawdata(:,2) = curcellrawdata(:,2)-curcellrawdata(1,2);
            % concatenating cell data to rawdata
            rawdata = [rawdata; curcellrawdata];
        end
        prevfolder = curfolder;
    end
    xlswrite(analysisplan,1,'experiments',['F' num2str(listnums(numofseries))]);
end

cd(output)
save([curfolder '_rawdata'],'rawdata');


end

end