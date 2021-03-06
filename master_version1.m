close all;
clearvars;
clc;

% Let user select analysis plan
cd '\\phys34212\migrationdata\MigrationData\Migration1\Repositories\cell-migration-analysis\analysisplan'
[file,path] = uigetfile('*.xls');
analysisplan = fullfile(path,file);
if isequal(file,0)
   disp('User selected Cancel');
end

[~,mypath,~] = xlsread(analysisplan,'path');
nd2path = mypath{1,2};
trackpath = mypath{2,2};
output = mypath{3,2};
bfmatlab = mypath{4,2};
addpath(bfmatlab);
featurespath = mypath{5,2};
codepath = mypath{6,2};
addpath(codepath);

% Extract metadata from nd2 files and save to output path
extractmetadata (nd2path, output, analysisplan);

% Extract raw data from tracking data
% exTRACKtRAWdata ( trackpath , output ); % uncomment to process in batches
% tic
% disp('Extracting raw data...')
extractRawData_alt (trackpath, output, analysisplan); % process each row individually
% disp('Raw data extracted!')
% toc

% Extract data from rawdata and metadata
% tic
% disp('Extracting data...')
extractData(output, analysisplan)
% disp('Data extracted!')
% toc

% Extract features from data
% tic
% disp('Extracting features...')
extractFeatures(output, analysisplan)
% disp('Features extracted!')
% toc

% Visualize data
% migrationvis(featurespath,output)
% tic
% disp('Visualizing features')
visualizeFeatures(output, analysisplan);
% disp('Features visualized!')
% tic

cd(codepath);