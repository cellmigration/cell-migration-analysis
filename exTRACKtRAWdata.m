function [  ] = exTRACKtRAWdata ( trackpath , output)
%Reads tracking data from mat file(s), extracts x,y location of each
%tracked cell in each stage (series). The output is a matrix where each row
%is defined as following: 1)relative y-coordinate 2)relative x-coordinate
%3)frame number 4)cell ID 5)series number

%   trackpath is the path where the mat files are found and outputpath is
%   where the raw variable will be saved. Each tracking data file should
%   be named based on this template:'yyyymmddCOL123C123S12' (it is the nd2
%   file name + 'S01' that represents the stage (series) number.

% list of experiments
analysisplan = 'D:\PhD Physics GT\Curtis Lab\Computation\Code files\analysisplan';
[num,~,raw] = xlsread(analysisplan,'experiments');
listnums = find(num(:,6)==0)+1;

% for seriesnum = 1:length(listnums)
%     curfolder = num2str(raw{listnums(seriesnum),1});
%     cd([trackpath '\' curfolder]);
    
    
for foldernum = 1:length(listnums)
    curfolder = num2str(raw{listnums(foldernum),1});
    cd([trackpath '\' curfolder]);
    filenames = dir;
    filenames = filenames(3:end-1);
    filenames = {filenames.name};
    rawdata = [];
    for numofseries = 1:size(filenames,2)
        curseries = filenames{numofseries};
        load(curseries);
        % adding series number in last column
        seriesnum = str2double(curseries(end-5:end-4));
        cellnums = unique(Mdati(:,4));
        Mdati = [Mdati, seriesnum*ones(size(Mdati,1),1)];
            for curcellnum = 1:length(cellnums)
                curcellind = find(Mdati(:,4)==cellnums(curcellnum));
                curcellrawdata = Mdati(curcellind,:);
                % changing coordinates relative to origin of cell
                curcellrawdata(:,1) = -curcellrawdata(:,1)+curcellrawdata(1,1);
                curcellrawdata(:,2) = curcellrawdata(:,2)-curcellrawdata(1,2);
                % concatenating cell data to rawdata
                rawdata = [rawdata; curcellrawdata];
            end
    end
    cd(output)
    save([curfolder '_rawdata'],'rawdata');
    xlswrite(analysisplan,1,'experiments',['F' num2str(listnums(foldernum))]);
end

if(isempty(listnums))
    fprintf('Tracking data for experiments in the list were already extracted\n')
end
end