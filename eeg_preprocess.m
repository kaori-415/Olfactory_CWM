function [Response] = eeg_preprocess(filename)

% preprocessing EEG data
% filtering, getting triggers, make epochs
% filtering: causal FIR bandpath filter and notch filter

%root = cd;
% addpath(cd)
% 
% [filename, path] = uigetfile('*.CSV');
% cd(path)

% definitions
Ch_order = {'Fz', 'Cz', 'Pz', 'Oz', 'F3', 'F4', 'P3', 'P4','A1', 'A2', 'EOG1', 'EOG2', 'trig'};
Fs = 1000; % sampling frequency 

event = [4352, 4353]; % event values
pretime = 500; posttime = 1000; %for ERP, .ms

filtorder = 20;
lpass = 50; hpass = 0.1;

condition = {'control', 'odorant'};

% data import
Data = importdata(filename);

% Reset Reference (A1A2)
A1 = Data.data(:,11);
A2 = Data.data(:,12);
ref = (A1+A2)/2;
Data.data(:,3:10) = Data.data(:,3:10) - repmat(ref,1,8);

timepoint = Data.data(:,1);
eventpoint = Data.data(:,2);

% filter
Data.filtdata = Data.data;

[Bf,a] = fir1(filtorder, [hpass/Fs, lpass/Fs]); % design a FIR filter
for nd = 3:size(Data.filtdata,2)-1
    Data.filtdata(:,nd) = filtfilt(Bf, a, Data.data(:, nd)); % FIR filter
    Data.filtdata(:,nd) = notch(Data.filtdata(:,nd)); % notch filter
end

% get event
trigger = Data.filtdata(:, end);
dtrigger = trigger(2:end)-trigger(1:end-1);
idx_d = find(dtrigger>15000)+1;
didx_d =idx_d(2:end)-idx_d(1:end-1);
idx_didx_d = find(didx_d < 100) + 1;
idx_d(idx_didx_d) = [];                %delete overlap trigger
idx = idx_d;
%length(idx)

inv_idx_d = find(dtrigger>15000);
inv_didx_d = inv_idx_d(2:end)-inv_idx_d(1:end-1);
inv_idx_didx_d = find(inv_didx_d <100);
inv_idx_d(inv_idx_didx_d) = [];
inv_idx = inv_idx_d;

% Seperate data by conditions
StartPoint = find(eventpoint == event(1));
EndPoint = find(eventpoint == event(2));


% check trigger
idx(min(StartPoint) > idx) = [];
idx(max(EndPoint) < idx) = [];
inv_idx(min(StartPoint) > inv_idx) = [];
inv_idx(max(EndPoint) < inv_idx) = []; 

% epoching data by memory onset
% time x chan x trial

memonset_idx = 1:3:length(idx);
memonset_trg = idx(memonset_idx);
n=1; ntrial = 1; ntrial_cond = 1;

while true
    
    if (memonset_trg(ntrial)) < StartPoint(n)
        ntrial = ntrial+1;
        continue
    end
    
    Memory(n).epoch(:,:,ntrial_cond) = ...
        Data.filtdata(memonset_trg(ntrial)-pretime:memonset_trg(ntrial)+posttime, [3:14]);
    
    ntrial = ntrial+1;
    ntrial_cond = ntrial_cond+1;
    
    if ntrial>length(memonset_trg)
        break
    end
    
    if memonset_trg(ntrial) >= EndPoint(n)
        n = n+1;
        ntrial_cond = 1;
    end
    
    if n > length(condition)
        break
    end
end

% epoching data by response onset
% time x chan x trial

for nresponse = 1:2
    
    response_idx = nresponse+1:3:length(idx);
    onset_trg = idx(response_idx);
    offset_trg = inv_idx(response_idx);
    
    n = 1;
    ntrial = 1;
    ntrial_cond = 1;
    
    while true
        
        Response(n).order(nresponse).epoch{ntrial_cond} = ...
            Data.filtdata(onset_trg(ntrial)-pretime:offset_trg(ntrial), [3:14]);
        Response(n).order(nresponse).rt(ntrial_cond) = (offset_trg(ntrial)-onset_trg(ntrial)+1)/Fs;
        
        ntrial = ntrial+1;
        ntrial_cond = ntrial_cond+1;
        
        if ntrial>length(onset_trg)
            break
        end
        
        if onset_trg(ntrial) > EndPoint(n)
            n = n+1;
            ntrial_cond = 1;
        end
        
        if n>length(condition)
            break
        end
    end
end

save PreEEG.mat Data Memory Response Ch_order

%cd(root)