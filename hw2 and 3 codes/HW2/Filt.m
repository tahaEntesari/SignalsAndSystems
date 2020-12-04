function [filtered] = Filt(x)
%
fc=40;
fs=120;
wc=2*pi*fc/fs;
n=-3:1/fs:3;
h=sin(n*wc)./(pi*n);
plot(n,h);
fftplot(h,1/fs,3,0,1);

end

