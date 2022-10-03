% Turning lvm files into txt files cause pythton
% doen't have a thing to do it

% danyon Stitt

clc
clear
close all

%% MULTIPLE FILE LVM TO TXT

myFolder = 'E:\bug fix octobver 2\';

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.lvm'); % Change to whatever pattern you need..
theFiles = dir(filePattern);

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Importing the file
    NAP = lvm_import(fullFileName);
    data = NAP.Segment1.data;
    filename_out = strrep(fullFileName, ".lvm", ".txt");
    writematrix(data, filename_out)
end

"complete"

%% SINGLE FILE LVM TO TXT

% full_file_name = "Steel no neck 30cm forehead no headgear.lvm";
% 
% NAP = lvm_import(full_file_name);
% data = NAP.Segment1.data;
% filename_out = strrep(full_file_name, ".lvm", ".txt");
% writematrix(data, filename_out)






