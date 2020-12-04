function [] = fftplot (x,Fs,time,which,newfig,xmin,xmax,ymin,ymax)
%   x: the vector
%   Fs:the sampling frequency
%   time:total lenght of the signal
%   which:  which=1   means that the given x is in frequency domain
%               which=0    means that you must work with fft(x)
%newfig:    whether or not to draw in a new figure.
  %             1 mean yes.0 no.

  
      
L=time*Fs;
if(which==1)
    P2 = abs(x/L);
elseif which==0
    P2=abs(fft(x)/L);
end
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
if(newfig==1)
figure;
end
plot(f,P1) ;
if nargin==9
    axis([xmin,xmax,ymin,ymax]);
end

end

