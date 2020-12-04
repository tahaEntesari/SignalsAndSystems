function [cumsumed_energy] = fftenergy(x,which,total)
% x is the vector of whose energy is to be calculated
% which=0 means x is in time domain;
%if total=1, the function will return the TOTAL energy rather than the
%whole power vector
if(which==0)
    x=fft(x);
end
cumsumed_energy=cumsum(x.*conj(x));
if(total==1)
    cumsumed_energy=cumsumed_energy(end);
end
end

