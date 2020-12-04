function [] = fftplot (input_vector ,Fs ,time ,which_domain ,newfig ,xmin ,xmax ,ymin ,ymax)
%   
%   Fs:the sampling frequency
%   time:total lenght of the signal
%   which_domain:  1   means that the given input is in frequency domain
%                               0  time domain
%   newfig:    whether or not to draw in a new figure.
%                   1 yes, 0 no.

  
      
L=time*Fs;
if(which_domain==1)
    P2 = abs(input_vector/L);
elseif which_domain==0
    P2=abs(fft(input_vector)/L);
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

