%delta 0.5-4Hz
y_delta1 = 2*4*sinc(4*t);
y_delta2 = 2*0.5*sinc(0.5*t);
y_delta = y_delta1/polyval(y_delta1,1) - y_delta2/polyval(y_delta2,1);
%%  filtering all signals delta region <4Hz
eeg_signal_filtered_delta = zeros(64,360,49);
for i=1:64
    for j=1:49
        eeg_signal_filtered_delta(i,:,j) = filter(y_delta,1,eeg_signal_ac_120(i,:,j));
    end
end
%%
for i=15 : 18
    figure
    SF = eeg_signal_filtered_delta(i,:,20);
    SF1 = eeg_signal_filtered_theta(i,:,20);
    SF2 = eeg_signal_filtered_alpha(i,:,20);
    SF3 = eeg_signal_filtered_beta(i,:,20);
    SF4 = eeg_signal_filtered_gamma(i,:,20);
Fs = 120;
Ts=1/Fs;
SFdft = fftshift(fft(SF));
df = Fs/length(SF);
freqvec = -Fs/2+df:df:Fs/2;
subplot(5,1,1)
plot(freqvec,2*abs(SFdft)/length(SF))
xlabel('Hz');
xlim([0, 40])
title('delta band')

SF1dft = fftshift(fft(SF1));
df = Fs/length(SF1);
freqvec = -Fs/2+df:df:Fs/2;
subplot(5,1,2)
plot(freqvec,2*abs(SF1dft)/length(SF1))
xlabel('Hz');
xlim([0, 40])
title('theta band')

SF2dft = fftshift(fft(SF2));
df = Fs/length(SF2);
freqvec = -Fs/2+df:df:Fs/2;
subplot(5,1,3)
plot(freqvec,2*abs(SF2dft)/length(SF2))
xlabel('Hz');
xlim([0, 40])
title('alpha band')

SF3dft = fftshift(fft(SF3));
df = Fs/length(SF3);
freqvec = -Fs/2+df:df:Fs/2;
subplot(5,1,4)
plot(freqvec,2*abs(SF3dft)/length(SF3))
xlabel('Hz');
xlim([0, 40])
title('beta band')

SF4dft = fftshift(fft(SF4));
df = Fs/length(SF4);
freqvec = -Fs/2+df:df:Fs/2;
subplot(5,1,5)
plot(freqvec,2*abs(SF4dft)/length(SF4))
xlabel('Hz');
xlim([0, 40])
title('gamma band')
end
%% test=zeros(16,4,360);
for i=1 : 16
    test(i,1,:)=eeg_signal_filtered_alpha(14+i,:,20);
        test(i,2,:)=eeg_signal_filtered_beta(14+i,:,20);
            test(i,3,:)=eeg_signal_filtered_gamma(14+i,:,20);
                test(i,4,:)=eeg_signal_filtered_delta(14+i,:,20);
end

save('test.mat', 'test');