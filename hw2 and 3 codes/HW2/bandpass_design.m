function [bb,ab] = bandpass_design(w1,w2)
% bandpass filter by truncating the impulse response
% w1 and w2 are the cutoff frequncys
% bb is the numerator
%           ************* w1 and w2 must be normalized*************
% there might be the need that the user finally normalize the gain by
% dividing the bb onto sth 


N=1400;
ah=zeros(1,floor(N/2));
ah(1)=1;
bh=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bh(i+1+floor(N/2))=exp(1j*(pi*(1/2+w1))*i)*sinc(i*1/2);
end
%bh=bh/3;
figure;
freqz(bh,ah)
title('highpass frequency responce')
%lowpass
N=200;
al=zeros(1,floor(N/2));
al(1)=1;
bl=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bl(i+1+floor(N/2))=sinc(i*w2);
end
bl=bl/3;
figure;
freqz(bl,al)
title('lowpass frequency response')
% bandpass

bb=conv(bl,bh);
ab=zeros(1,length(bb));
ab(1)=1;
figure
freqz(bb,ab)
title('bandpass frequency response');
end

