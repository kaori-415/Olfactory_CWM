function y = notch(x)
% Notch Filtering 

%
% MATLAB Code
% Developed by MATLAB(R) 7.14 & DSP System Toolbox 8.2 
%
% 10-Sep-2012 17:43:20
%

persistent Hd;

if isempty(Hd)
    
    N  = 2;     % Order
    F0 = 60;    % Center frequency
    BW = 1.5;     % Bandwidth
    Fs = 1000;  % Sampling Frequency
    
    h = fdesign.notch('N,F0,BW', N, F0, BW, Fs);
    
    Hd = design(h, 'butter', ...
        'SOSScaleNorm', 'Linf');
    
    
    
    set(Hd,'PersistentMemory',true);
    
end

y = filter(Hd,x);


% [EOF]
