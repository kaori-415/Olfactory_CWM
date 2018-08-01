function Conditions = define_condition(index_control,condition,Response)

% define conditions by output from FourColorExp_hue.m

% get Response
if nargin < 3
    % check dir
    ex = exist('result.mat', 'file');
    
    if ex == 0
        dpath = uigetdir('Data/s01', 'Select folder including a result file');
        cd(dpath)
    end
    

end

if ~exist('Response', 'var')
    
    try
        Response = importdata('result.mat');
    catch
        error('cannot find a result.mat file');
    end
end

tmpresult = load('result.mat');
save tmpresult.mat tmpresult

% make a condition label

Conditions.labels = cell(1,size(Response,1));
Conditions.index = zeros(1,size(Response,1));

Conditions.labels{index_control} = 'control';
Conditions.index(index_control) = 1; 

index_conditions = 1:length(Conditions.labels);

if nargin < 2 || isempty(condition)
    
    odor_label = {'odorant'};
    if size(Response,1) > 2
        num_cell = num2cell(1:size(Response,1)-1);
        num_cell = cellfun(@num2str, num_cell, 'UniformOutput', false);
        
        odor_label = repmat(odor_label, size(num_cell));
        odor_label = strcat(odor_label, num_cell);
    end
    
else
    odor_label = condition;
    
    if length(condition) ~= length(Conditions)-1
        error('condition input is not matched')
    end
end
    
Conditions.labels(index_conditions ~= index_control) = odor_label(:);
Conditions.index(index_conditions ~= index_control) = 2:size(Response,1);

save result.mat Response Conditions
disp(Conditions)
