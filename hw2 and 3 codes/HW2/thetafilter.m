function y = thetafilter(x)
%DOFILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.3 and DSP System Toolbox 9.5.
% Generated on: 23-Mar-2018 20:11:28

persistent Hd;

if isempty(Hd)
    
    N      = 2000;  % Order
    Fstop1 = 3.5;   % First Stopband Frequency
    Fpass1 = 4;     % First Passband Frequency
    Fpass2 = 7;     % Second Passband Frequency
    Fstop2 = 7.5;   % Second Stopband Frequency
    Fs     = 120;   % Sampling Frequency
    
    h = fdesign.bandpass('n,fst1,fp1,fp2,fst2', N, Fstop1, Fpass1, Fpass2, ...
        Fstop2, Fs);
    
    Hd = design(h, 'equiripple');
    
    
    
    set(Hd,'PersistentMemory',true);
    
end

y = filter(Hd,x);

