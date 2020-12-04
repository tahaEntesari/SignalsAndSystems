
%% Question 1
clear

load G_N_I.mat
Fs=2*10^5;
number=numel(note);
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)

%audiowrite('music.wav',m1,Fs);
%% dc addition to note
note=note+10;
Fs=2*10^5;
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)
note=note-10;
%% dc addition to interval
interval=interval+.1;
Fs=2*10^5;
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)
interval=interval-.1;


%% noise to note
load G_N_I.mat
Fs=2*10^5;
note=note+randn(1,number);
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)
%% noise to interval
load G_N_I.mat
Fs=2*10^5;
interval=interval+randn(1,number);
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)
%% noise to gain
load G_N_I.mat
Fs=2*10^5;
gain=gain+randn(1,number);
m1=SoundMaker(gain,note,interval,Fs);
sound(m1,Fs)
%% noise to m1 with variance 1
load G_N_I.mat
Fs=2*10^5;
m1=SoundMaker(gain,note,interval,Fs);
m1=m1+randn(1,numel(m1));
sound(m1,Fs)
%% maximum power noise
load G_N_I.mat
Fs=2*10^5;
m1=SoundMaker(gain,note,interval,Fs);
m1=m1+20*randn(1,numel(m1));
sound(m1,Fs)
%% m1+m2 and there ffts
load G_N_I.mat
Fs=2*10^5;
m1=SoundMaker(gain,note,interval,Fs);
m2=m1+4*randn(1,numel(m1));
   sound(m2,Fs)
%                   audiowrite('m2 music.wav',m2,Fs);
%
time=sum(interval);
L=time*Fs;
m1ft=fft(m1);
m2ft=fft(m2);

spec1nonoise=abs(m1ft/(time*Fs));
spec2nonoise=spec1nonoise(1:L/2+1);
spec2nonoise(2:end-1)=2*spec2nonoise(2:end-1);

spec1noise=abs(m2ft/(time*Fs));
spec2noise=spec1noise(1:L/2+1);
spec2noise(2:end-1)=2*spec2noise(2:end-1);

f=Fs*(0:(L/2))/L;

figure(1);
subplot(1,2,1),plot(f,spec2nonoise),title('The main signal with no noise'),axis([100,500,0,.45]),xlabel('frequency'),ylabel('gain');
subplot(1,2,2),plot(f,spec2noise),title('Signal with noise with STD 2'),axis([100,500,0,.45]),xlabel('frequency'),ylabel('gain');
figure(2);
plot(f,abs(spec2nonoise-spec2noise)),title('absolute error of the two signals'),axis([0,1000,0,.01]),xlabel('frequency'),ylabel('gain difference');

%filtering

%% sinc conv
spmd
t=0:1/Fs:(sum(interval));
sinced=sinc(400*t);
tic
m2_f=conv(sinced,m2);
toc
m2_f=m2_f(1:numel(m2))/400;
sound(m2_f,Fs)
end

audiowrite('m2_f.wav',m2_f,Fs);

%% trying to reach better noise

%% exp*sinc 
load G_N_I.mat
Fs=2*10^5;
m1=SoundMaker(gain(1:10),note(1:10),interval(1:10),Fs);
m2=m1+4*randn(1,numel(m1));
sinced=2*cos(200*t).*sinc(100*t);
tic
m2_best=conv(sinced,m2);
toc
%%
sinced=sinc(400*t)/400-sinc(200*t)/200;
tic
m2_2=conv(sinced,m2);
toc
sound(real(m2_2),Fs)
%%
sinced=exp(200j*t).*sinc(200*t);
tic
m2_3=conv(sinced,m2);
toc
sound(real(m2_3),Fs)

%%
%%
%%
%% question 2 part a
n=1000;
w=-10:n;
x=zeros(1,n);
x(11:16)=1;
h=exp(-w).*(heaviside(w)-heaviside(w-5));
figure(3);
cla;
stem(conv(x,h));
ylabel('convolution');
hold on
stem(MyConv(x,h));
axis([10,50,0,1.5])
legend('matlab conv','my conv');
%% part b
t=-10:1/20:20;
x=-heaviside(t+1)+3*heaviside(t)-2*heaviside(t-2);
h=heaviside(t)-heaviside(t-1).*heaviside(t)-heaviside(t-1);
figure(4);
cla;
stem(conv(x,h));
ylabel('convolution');
hold on
stem(MyConv(x,h));
legend('matlab conv','my conv');
axis([300,900,-80,60])
%% part c
x=sin(2*pi*t)+sin(2*pi*100*t);
h=sinc(2*t);
h(401:end)=0;
figure(5);
cla;
stem(conv(x,h));
ylabel('convolution');
hold on
stem(MyConv(x,h));
legend('matlab conv','my conv');
axis([0,1000,-10,10])
%% part d
x=zeros(1,numel(w));
h=zeros(1,numel(w));
x(7)=-1;
x(11)=1;
x(15)=2;
h(11)=1;
figure(6);
cla;
stem(conv(x,h));
ylabel('convolution');
hold on
stem(MyConv(x,h));
legend('matlab conv','my conv');
axis([0,100,-1,2])


