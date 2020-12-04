%%  Question 1

    global Fs time L;
    Fs=2400;
    time=3;
    L=Fs*time;
    f = Fs*(0:(L/2))/L;
    load EEG_signal
    [row,column,dep]=size(EEG_signal_edited);
    for i=1:5
            fftplot(EEG_signal_edited(1,:,i),Fs,time,0,1,0,60,-inf,inf)
            txt=sprintf('frequency domain  activity %d for the first channel',i);
            title(txt);
            xlabel('frequency');
            ylabel('gain');
    end
    %% finding cutoff freq
    totenergy=zeros(5,7200);
    for i=1:5
        ffted(i).fft=fft(EEG_signal_edited(1,:,i));
    end
    index_90percent=zeros(1,5);
    for i=1:5
        totenergy(i,:)=cumsum(abs(ffted(i).fft.^2));
        a=find(totenergy(i,:)>.9*totenergy(i,end));
        index_90percent(i)=a(1);
    end
    clear ffted

    
    %% question 1 edited by deleting mean
    b=mean(EEG_signal_edited(:,:,:),2);
    EEG_signal=EEG_signal_edited-b;
    clear b;
    
     for i=1:5
            fftplot(EEG_signal(1,:,i),Fs,time,0,1,0,60,-inf,inf)
            txt=sprintf('frequency domain  activity %d for the first channel',i);
            title(txt);
            xlabel('frequency');
            ylabel('gain');
    end
    %% finding cutoff freq after canceling signal mean
    totenergy=zeros(5,3601);
    for i=1:5
        ffted(i).fft=fft(EEG_signal(1,:,i));
        P1 = ffted(i).fft(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
    end
    index_90percent=zeros(1,5);
    for i=1:5
        totenergy(i,:)=cumsum(abs(P1.^2));
        a=find(totenergy(i,:)>.9*totenergy(i,end));
        index_90percent(i)=a(1);
    end
clear ffted
%%  down sampling
Fs_down=120;
L_down=Fs_down*time;
down_EEG=EEG_signal(:,1:20:end,:);
f_down=Fs_down*(0:(L_down/2))/L_down;
%% cutoff freq after downsapmling
totenergy=zeros(5,181);
    for i=1:5
        ffted(i).fft=fft(down_EEG(1,:,i));
        P1 = ffted(i).fft(1:L_down/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
       
    end
    index_90percent=zeros(1,5);
    for i=1:5
        totenergy(i,:)=cumsum(abs(P1.^2));
        a=find(totenergy(i,:)>.9*totenergy(i,end));
        index_90percent(i)=a(1);
    end
clear ffted


%% design filter by zeros
r=[.25,.5,.9,.95,1.05,1.1,3];
w=pi/2;
N=17;
a=zeros(1,N);
a(1)=1;
Z=zeros(1,N);
syms z
for j=1:7
 figure
 for i=1:N
     Z(i)=r(j)*exp(1j*(w+abs(2*pi-2*w)*i/N)); 
 end
 H=@(z)prod(z-Z);
sth=expand(H(z));
b=eval(coeffs(sth));
freqz(b,a);
atxt=sprintf('r=%f, w=%f,N=%d',r(j),w,N);
title(atxt)
end




%% filter design by truncating impulse response

a=zeros(1,5);
a(1)=1;
b=zeros(1,8);
for i=-3:4
    b(i+4)=sinc(i/2);
end
b=b/2.5;
freqz(b,a)

%% filtering
w=35/120;
N=27;
a=zeros(1,floor(N/2));
a(1)=1;
b=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    b(i+1+floor(N/2))=sinc(i*w);
end
b=b/3;
figure;
freqz(b,a)
txt=sprintf('filter constructed by N=%d,truncating impulse response',N);
title(txt);
figure
[b_butter,a_butter]=butter(N,w);
freqz(b_butter,a_butter);
txt=sprintf('filter constructed by N=%d,butterworth',N);
title(txt);
filteredEEG_trunc=filter(b,a,down_EEG(64,:,5));
filteredEEG_butter=filter(b_butter,a_butter,down_EEG(1,:,1));

 P2_trunc=abs(fft(filteredEEG_trunc)/L_down);
 P2_butter=abs(fft(filteredEEG_butter)/L_down);
P1_trunc = P2_trunc(1:L_down/2+1);
P1_trunc(2:end-1) = 2*P1_trunc(2:end-1);
energy_trunc=cumsum(P1_trunc.^2);

P2_butter=abs(fft(filteredEEG_butter)/L_down);
P1_butter= P2_butter(1:L_down/2+1);
P1_butter(2:end-1) = 2*P1_butter(2:end-1);
energy_butter=cumsum(P1_butter.^2);


a=find(energy_trunc>.999*energy_trunc(end));
trunced_index=a(1);
a=find(energy_butter>.999*energy_butter(end));
butter_index=a(1);
fprintf('frequency for which 99.9 percent of the energy\nlies before it is:\n 1.truncated impulse response %f\n2.butterworth filter %f\n',f_down(trunced_index),f_down(butter_index))
%%
fftplot(filteredEEG_trunc,120,3,0,1)
title('frequency domain of filtered signal by truncated impulse response');

fftplot(filteredEEG_butter,120,3,0,1)
title('frequency domain of filtered signal by butterworth');
figure
subplot(1,2,1)
plot(down_EEG(1,:,1));
title('unfiltered');
subplot(1,2,2)
plot(filteredEEG_trunc);
title('truncation filtering result');
figure
subplot(1,2,1)
plot(down_EEG(1,:,1));
title('unfiltered');
subplot(1,2,2)
plot(filteredEEG_butter);
title('butter filtering result');

%% highpass design
w=1/4;
N=48;
a=zeros(1,floor(N/2));
a(1)=1;
b=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    b(i+1+floor(N/2))=exp(1j*1*pi*i)*sinc(i*w);
end
b=b/3;
figure;
freqz(b,a)
title('high pass filter by truncating impulse response');
%% bandpass
%highpass
w=1/2;
N=20000;
ah=zeros(1,floor(N/2));
ah(1)=1;
bh=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bh(i+1+floor(N/2))=exp(1j*(pi*(1/2+35/120))*i)*sinc(i*w);
end
bh=bh/3;
figure;
freqz(bh,ah)
title('hp')
%lowpass
w=1/3;
N=20000;
al=zeros(1,floor(N/2));
al(1)=1;
bl=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bl(i+1+floor(N/2))=sinc(i*w);
end
bl=bl/3;
figure;
freqz(bl,al)

bb=conv(bl,bh);

ab=zeros(1,length(bb));
ab(1)=1;
figure
freqz(bb,ab)






    %% Question 2
    %%  part 1
    for i=1:3
        activity(i).mean=mean(EEG_signal_edited(:,:,i),2);
        activity(i).std=std(EEG_signal_edited(:,:,i),0,2);
        activity(i).max=max(EEG_signal_edited(:,:,i),[],2);
        activity(i).min=min(EEG_signal_edited(:,:,i),[],2);
    end

    %%  part 2
    t=linspace(0,time,7200);
        for i=1:5
            for j=1:49
                a=sprintf(' chanel %d  activity %d',i,j);
                figure;
                subplot(1,2,1);
                plot(t,EEG_signal(i,:,j))
                title(a);
                subplot(1,2,2);
                fftplot(EEG_signal(i,:,j),Fs,time,0,0,0,70,-inf,inf);
                
            end
        end
%% part 4 filling outliers
for i=1:64
    for j=1:49
        down_EEG_outlied(i,:,j)=filloutliers(down_EEG(i,:,j),'spline','gesd');
    end
end
 %%     filtering the signal
 w=35/120;
N=61;
a=zeros(1,floor(N/2));
a(1)=1;
b=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    b(i+1+floor(N/2))=sinc(i*w);
end
b=b/3;
N=27;
[b_butter,a_butter]=butter(N,w);
filteredEEG_trunc=zeros(64,360,49);
filteredEEG_butter=zeros(64,360,49);
for i=1:64
    for j=1:49
filteredEEG_trunc(i,:,j)=filter(b,a,down_EEG_outlied(i,:,j));
filteredEEG_butter(i,:,j)=filter(b_butter,a_butter,down_EEG_outlied(i,:,j));
    end
end


%%

for j=1:49
    figure
subplot(2,2,1)
plot(filteredEEG_trunc(1,:,j));
title('truncation filtering');
subplot(2,2,2)
plot(filteredEEG_butter(1,:,j));
title('butterworth filtering');
subplot(2,2,[3,4]);
plot(down_EEG_outlied(1,:,j));
txt=sprintf('the unfiltered signal channel 1 activity %d',j);
title(txt);

end
        %% part 5
        % since no particular pattern could be found for the numbers of 
        % corresponding left and right channels,I have entered them here by hand!!!!
        common={[2,3],[34,35],[4,5],[6,10],[36,39],[7,9],[37,38],[40,44],[11,14],[41,43],[12,13],[15,19],...
            [45,48],[16,18],[46,47],[49,53],[20,23],[50,52],[21,22],[24,28],[54,57],[25,27],[55,56],[29,30],...
            [58,60],[31,32],[62,63],[64,33]};
        %choose n for the number of the Activity;
        n=27;
        
        for j=1:28
                figure
                subplot(2,2,1);
                plot(EEG_signal_edited(common{j}(1),:,n),EEG_signal_edited(common{j}(2),:,n));
                txt=sprintf('plot of channel %d with respect to channel %d',common{j}(2),common{j}(1));
                title(txt);
                txt=sprintf('channel %d',common{j}(1));
                xlabel(txt);
                txt=sprintf('channel %d',common{j}(2));
                ylabel(txt)
                subplot(2,2,3);
                plot(abs(EEG_signal_edited(common{j}(1),:,n)-EEG_signal_edited(common{j}(2),:,n)));
                title('Absolute difference plot');
                subplot(2,2,2)
                txt=sprintf('activity %d plots for channels %d and %d',n,common{j}(1),common{j}(2));
                text(0,0,txt,'interpreter','latex'); 
                axis off;
               
            end
        
        %% part 6
        [correlation,pval]=corr(EEG_signal_edited(:,:,1)');
        [x_index,y_index]=find(correlation >.95 & correlation<1);
        
        %% part 7
        figure;
        boxplot(EEG_signal_edited(1,:,1)')
        title('boxplot of first channel');
        figure;
        boxplot(EEG_signal_edited(1:16,:,1)')
        title('boxplot of channels 1 to 16 for the first activity');

        
        %% part 7 with edited mean
        figure;
        boxplot(EEG_signal(1,:,1)')
        title('boxplot of first channel with zero mean');
        figure;
        boxplot(EEG_signal(1:16,:,1)')
        title('boxplot of channels 1 to 16 for the first activity with zero mean');
        figure;
        boxplot(filteredEEG_trunc(1:16,:,1)');
        title('boxplot of channels 1 to 16 for the first activity with filtered signal');
        %% part 8
    % alpha band is 8-15
    [b_alpha,a_alpha]=bandpass_design(8/120,15/120);
    b_alpha=180*b_alpha/6;
    title(' alpha bandpass filter');
    % beta is 14 til 31
    [b_beta,a_beta]=bandpass_design(14/120,31/120);
    b_beta=10*b_beta/3;
    
    title(' beta bandpass filter');
    % gamma is anything bigger than 32
    [b_gamma,a_gamma]=bandpass_design(32/120,1);
    b_gamma=10*b_gamma*1.5;
    title(' gamma bandpass filter');
    %delta is anything lower than 4 
    [b_delta,a_delta]=bandpass_design(0,6/120);
     b_delta=100*b_delta/12;
    title(' delta bandpass filter');
%close all
    for i=1:16
        
        band(i).alpha=filter(b_alpha,a_alpha,(filteredEEG_butter(i,:,1)));
        band(i).beta=filter(b_beta,a_beta,filteredEEG_butter(i,:,1));
        band(i).gamma=filter(b_gamma,a_gamma,filteredEEG_butter(i,:,1));
        band(i).delta=filter(b_delta,a_delta,filteredEEG_butter(i,:,1));
        
    end
    %% filter by matlab
        for i=1:16
        
        band_m(i).alpha=alphafilter((filteredEEG_butter(i,:,1)));
        band_m(i).beta=betafilter(filteredEEG_butter(i,:,1));
        band_m(i).gamma=gammafilter(filteredEEG_butter(i,:,1));
        band_m(i).delta=deltafilter(filteredEEG_butter(i,:,1));
        
    end
    
    %% plots using my own filter
    channels =[2,5,7,16];
    for i=channels
    figure
    subplot(2,5,1);
    fftplot(filteredEEG_butter(i,:,1),120,3,0,0);
    txt=sprintf('frequency domain of normal signal \nchannel %d',i);
    title(txt);
    
    subplot(2,5,2);
    fftplot(band(i).alpha,120,3,0,0);
    txt=sprintf('frequency domain of alpha band');
    title(txt);
    
    subplot(2,5,6)
    plot(filteredEEG_butter(i,:,1));
    txt=sprintf('time domain of normal signal ');
    title(txt);
    
    subplot(2,5,7)
    plot(1:360,(band(i).alpha));
    txt=sprintf('time domain of alphaband signal');
    title(txt);
    
    subplot(2,5,3);
    fftplot(band(i).beta,120,3,0,0);
    txt=sprintf('frequency domain of beta band');
    title(txt);
    
    subplot(2,5,8)
    plot(1:360,band(i).beta);
    txt=sprintf('time domain of beta band signal');
    title(txt);
    
    subplot(2,5,4);
    fftplot(band(i).gamma,120,3,0,0);
    txt=sprintf('frequency domain of gamma band');
    title(txt);
    
    subplot(2,5,9)
    plot(1:360,band(i).gamma);
    txt=sprintf('time domain of gamma band signal');
    title(txt);
    
    subplot(2,5,5);
    fftplot(band(i).delta,120,3,0,0);
    txt=sprintf('frequency domain of delta band');
    title(txt);
    
    subplot(2,5,10)
    plot(1:360,band(i).delta);
    txt=sprintf('time domain of delta band signal');
    title(txt);
    end
    %% plots using matlab filter
    for i=channels
    figure
    subplot(2,5,1);
    fftplot(filteredEEG_butter(i,:,1),120,3,0,0);
    txt=sprintf('frequency domain of normal signal \nchannel %d',i);
    title(txt);
    
    subplot(2,5,2);
    fftplot(band_m(i).alpha,120,3,0,0);
    txt=sprintf('frequency domain of alpha band');
    title(txt);
    
    subplot(2,5,6)
    plot(filteredEEG_butter(i,:,1));
    txt=sprintf('time domain of normal signal ');
    title(txt);
    
    subplot(2,5,7)
    plot(1:360,(band_m(i).alpha));
    txt=sprintf('time domain of alphaband signal');
    title(txt);
    
    subplot(2,5,3);
    fftplot(band_m(i).beta,120,3,0,0);
    txt=sprintf('frequency domain of beta band');
    title(txt);
    
    subplot(2,5,8)
    plot(1:360,band_m(i).beta);
    txt=sprintf('time domain of beta band signal');
    title(txt);
    
    subplot(2,5,4);
    fftplot(band_m(i).gamma,120,3,0,0);
    txt=sprintf('frequency domain of gamma band');
    title(txt);
    
    subplot(2,5,9)
    plot(1:360,band_m(i).gamma);
    txt=sprintf('time domain of gamma band signal');
    title(txt);
    
    subplot(2,5,5);
    fftplot(band_m(i).delta,120,3,0,0);
    txt=sprintf('frequency domain of delta band');
    title(txt);
    
    subplot(2,5,10)
    plot(1:360,band_m(i).delta);
    txt=sprintf('time domain of delta band signal');
    title(txt);
    end
    

    
    %% part 9       the required channels [2,5,7,16]
    % 4 channels
    %12 windows of 30 datas
    windowed_EEG=zeros(4,30,12);
    alpha=zeros(4,30,12);
    beta=zeros(4,30,12);
    gamma=zeros(4,30,12);
    delta=zeros(4,30,12);
    alphaenergy=zeros(4,12);
    betaenergy=zeros(4,12);
    gammaenergy=zeros(4,12);
    deltaenergy=zeros(4,12);
    
    k=1;
    % bandpass
%highpass
w=1/2;
N=1400;
ah=zeros(1,floor(N/2));
ah(1)=1;
bh=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bh(i+1+floor(N/2))=exp(1j*(pi*(1/2+8/120))*i)*sinc(i*w);
end
bh=bh/3;

%lowpass
w=15/120;
N=200;
al=zeros(1,floor(N/2));
al(1)=1;
bl=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bl(i+1+floor(N/2))=sinc(i*w);
end
bl=30*bl;
bb=conv(bl,bh);

ab=zeros(1,length(bb));
ab(1)=1;
    
    for i=channels
        for j=1:12
        windowed_EEG(i,:,j)=filteredEEG_butter(i,(30*(j-1)+1):(30*j),1);
        delta(k,:,j)=deltafilter(windowed_EEG(i,:,j));
        gamma(k,:,j)=gammafilter(windowed_EEG(i,:,j));
        beta(k,:,j)=betafilter(windowed_EEG(i,:,j));
        alpha(k,:,j)=filter(bb,ab,windowed_EEG(i,:,j));
   
        alphaenergy(k,j)=fftenergy(alpha(k,:,j),0,1);
        betaenergy(k,j)=fftenergy(beta(k,:,j),0,1);
        gammaenergy(k,j)=fftenergy(gamma(k,:,j),0,1);
        deltaenergy(k,j)=fftenergy(delta(k,:,j),0,1);
    
        end
        k=k+1;
    end
    for i=1:4
        figure
        subplot(2,2,1)
        stem(1:12,alphaenergy(i,:));
        txt=sprintf('alpha band energy for 250ms consequent windows\nchannel %d',channels(i));
        title(txt)
        
        subplot(2,2,2)
        stem(1:12,betaenergy(i,:));
        title('beta band energy for 250ms consequent windows');
        
        subplot(2,2,3)
        stem(1:12,gammaenergy(i,:));
        title('gamma band energy for 250ms consequent windows');
        
        subplot(2,2,4)
        stem(1:12,deltaenergy(i,:));
        title('delta band energy for 250ms consequent windows');
    end
    
        
    
    
%%  part 11 
% bandpass
%highpass
w=1/2;
N=1400;
ah=zeros(1,floor(N/2));
ah(1)=1;
bh=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bh(i+1+floor(N/2))=exp(1j*(pi*(1/2+8/120))*i)*sinc(i*w);
end
bh=bh/3;

%lowpass
w=15/120;
N=200;
al=zeros(1,floor(N/2));
al(1)=1;
bl=zeros(1,N);
for i=-floor(N/2):floor(N/2)
    bl(i+1+floor(N/2))=sinc(i*w);
end
bl=bl*30;
bb=conv(bl,bh);

ab=zeros(1,length(bb));
ab(1)=1;
alphaenergy=zeros(1,49);
gammaenergy=zeros(1,49);
betaenergy=zeros(1,49);
deltaenergy=zeros(1,49);
thetaenergy=zeros(1,49);
n=27;
     for i=1:49
        
        band_m(i).alpha=filter(bb,ab,(filteredEEG_butter(n,:,i)));
        band_m(i).beta=betafilter(filteredEEG_butter(n,:,i));
        band_m(i).gamma=gammafilter(filteredEEG_butter(n,:,i));
        band_m(i).delta=deltafilter(filteredEEG_butter(n,:,i));
        band_m(i).theta=thetafilter(filteredEEG_butter(n,:,i));
        
        alphaenergy(i)=sum(abs(band_m(i).alpha.^2));
        gammaenergy(i)=sum(abs(band_m(i).gamma.^2));
        betaenergy(i)=sum(abs(band_m(i).beta.^2));
        deltaenergy(i)=sum(abs(band_m(i).delta.^2));
        thetaenergy(i)=sum(abs(band_m(i).theta.^2));
        
     end
    energy=zeros(1,49);
    energy(:)=sum(filteredEEG_butter(n,:,:).^2);
    %%
    figure
    stem(alphaenergy)
    title('alpha');
    figure
    stem(gammaenergy)
    title('gamma');
    figure
    stem(betaenergy)
    title('beta');
    figure
    stem(deltaenergy)
    title('delta');
    %%
    figure
    subplot(2,1,1)
    stem(energy(1:49))
    title('channel 1 total energy');
    subplot(2,1,2)
    stem(thetaenergy)
    title('theta band energy');
    findex=find(thetaenergy>mean(thetaenergy));
    sindex=find(energy<mean(energy));

    