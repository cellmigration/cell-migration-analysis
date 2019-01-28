function [  ] = extractmetadata (nd2path, output, analysisplan)
%Reads metadata from nd2 file(s) and saves a metadata struct variable for
%each nd2 file. The saved variable is named as the file + '_metadata'
%
%   [ filenames ] = extractmetadata ( nd2path , outputpath )
%
%   nd2path is the path where the nd2 files are found and outputpath is
%   where the metadata variables will be saved.
%
%   metadata is a struct variable where metadata.um2px is micrometer to
%   pixel ratio and metadata.time is a matrix of timestamps in minutes.
%   Rows correspond to stages (series) and columns to timepoints.
%
%   This function uses "bfopen", hence it should be installed on the
%   computer and its path should be added.

%% read file names
[num,~,raw] = xlsread(analysisplan,'experiments');
%Find the excel index of series whose metadata has not been read yet
filenums = find(num(:,4)==0)+1;

% load analysisplan1
% filenums = find(analysisplan1(:,7)==0);

% read files
for curfilenum = 1:size(filenums,1)
    % read name of nd2 file; convert to string if necessary
    curfile = num2str(raw{filenums(curfilenum),2});
    curfile_out = num2str(raw{filenums(curfilenum),1});
    volume = bfopen(curfile);%bfOpen3DVolume(curfile);
    metadata.um2px = str2double(volume{1,2}.get('Global dCalibration'));
    numofseries = size(volume,1);
    numoftimepoints = volume{1,2}.get('Global Time Loop');
    if size(numoftimepoints,1)==0
        numoftimepoints = 1;
    else
        numoftimepoints = str2double(numoftimepoints);
    end
    % read series and timepoints of current file
    timetable = zeros(numofseries*numoftimepoints,1);
    ind = 1;
    for curseries = 1:numofseries
        for curtimepoint = 1:numoftimepoints
            if curtimepoint < 10
                if numoftimepoints > 9
                    curtimestr = ['timestamp #0',num2str(curtimepoint)];
                else
                    curtimestr = ['timestamp #',num2str(curtimepoint)];
                end
            else
                curtimestr = ['timestamp #',num2str(curtimepoint)];
            end
            curtime = volume{curseries,2}.get(curtimestr);
            timetable(ind) = curtime;
            ind = ind+1;
        end
    end
    savestr = [curfile_out,'_metadata'];
    timetable = reshape(timetable,numofseries,numoftimepoints)/60;
    metadata.time = timetable; % each row is a different series and columns are timepoints
    cd(output)
    save(savestr,'metadata')
    xlswrite(analysisplan,1,'experiments',['G' num2str(filenums(curfilenum))]);
end
if(isempty(filenums))
    fprintf('metadata for experiments in the list were already extracted\n')
end
end