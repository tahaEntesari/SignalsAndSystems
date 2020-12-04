w=1/30;
N=2000;
al=zeros(1,floor(N/2));
al(1)=1;
bl=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bl(i+1+floor(N/2))=sinc(i*w);
end
bl=bl/3;
figure;
freqz(bl,al)