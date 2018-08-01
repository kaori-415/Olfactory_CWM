function eeg_preprocess_auto(diridx)

root =cd; 
addpath(root);

resultsRoot = 'Data';
DirHeadStr = 's';

% get data
cd(resultsRoot);

% make dirnames
diridxCell = num2cell(diridx);
diridxStr = cellfun(@num2str, diridxCell, 'UniformOutput', false);
zeropadIdx = find(cellfun(@length, diridxStr)<2);
diridxStr(zeropadIdx) = strcat('0', diridxStr(zeropadIdx)); % padding zero as character like '01'
dirnames = strcat(DirHeadStr, diridxStr); % make dirnames including results
dirnamesLen = length(dirnames);

dirnames_withpath = fullfile(dirnames, 'EEG_source'); 

for n = 1:dirnamesLen
    
    cd(dirnames_withpath{n});
    fname = dir('*CSV');
    
    eeg_preprocess(fname.name);
    % eeg_preprocess_checktrig
    
    cd('../../');
end
    