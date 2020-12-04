function y = thetafilter(x)
%DOFILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.3 and DSP System Toolbox 9.5.
% Generated on: 03-Apr-2018 15:57:30

persistent Hd;

if isempty(Hd)
    
    Fstop1 = 3.4;  % First Stopband Frequency
    Fpass1 = 3.5;  % First Passband Frequency
    Fpass2 = 7.5;  % Second Passband Frequency
    Fstop2 = 7.6;  % Second Stopband Frequency
    Astop1 = 60;   % First Stopband Attenuation (dB)
    Apass  = 1;    % Passband Ripple (dB)
    Astop2 = 60;   % Second Stopband Attenuation (dB)
    Fs     = 120;  % Sampling Frequency
    
    h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
        Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
    
    Hd = design(h, 'equiripple', ...
        'MinOrder', 'any');
    
    
    
    set(Hd,'PersistentMemory',true);
    
end

y = filter(Hd,x);

