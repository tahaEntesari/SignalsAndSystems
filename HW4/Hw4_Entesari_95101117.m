
%% Question 1
clear
b1=[1,1];
a1=[1,2];
b2=[1];
a2=[500,0,0];
g1=tf(b1,a1);
g2=tf(b2,a2);
G=series(g1,g2);
fback=tf(1,1);
%1
H=feedback(G,fback);
%2
Poles=pole(H);
%3
figure;
pzmap(H)
%4
stable=isstable(H);
%5
figure;
impulse(H);
%6
figure;
bode(H);
%8
syms t 
x=@(t)exp(-t/10000).*cos(pi/500*t).*(t>0);
t=-5000:.1:5000;
X=x(t);
resp=lsim(H,X',t);
figure;
plot(t,resp);
xlabel('time');
title('response of the system, to the given input');
axis([-100,t(end),-inf,inf])

%% Question 2 part 1
clear
wn=1000;
u=[3,1.5,1,.5,0,-.5,-1.5];
for i=u
H=tf(wn^2,[1,2*i*wn,wn^2]);
figure
subplot(1,2,1);
bode(H);
txt=sprintf('Bode diagram, u=%f, wn=%d',i,wn);
title(txt);
subplot(1,2,2)
pzmap(H);
end
%% part 4
wn=100;
u=.7;
H=tf(wn^2,[1,2*u*wn,0]);
subplot(1,2,1);
bode(H)
txt=sprintf('Bode diagram, u=%f, wn=%d',u,wn);
title(txt);
subplot(1,2,2);
pzmap(H);
%% part 5
Tp=logspace(-5,100,20);
Tp=[0,Tp];
for i=Tp
    H=tf(1,[i,1]);
    figure
    subplot(2,3,[1,4]);
    bode(H);
    txt=sprintf('Bode diagram, Tp=%e',i);
    title(txt);
    subplot(2,3,[2,3]);
    pzmap(H);
    subplot(2,3,[5,6]);
    impulse(H);
end
%% part 7
    wn=100;
    u=.3;
    Tp=220;
    h2=tf(wn^2,[1,2*u*wn,0]);
    h1=tf(1,[Tp,1]);
    H=(series(h1,h2));
    figure
    subplot(2,3,[1,4]);
    bode(H);
    title('system without feedback');
    subplot(2,3,[2,3]);
    pzmap(H);
    subplot(2,3,[5,6]);
    impulse(H);
    %%
    fback=tf(1,1);
    H=feedback(series(h1,h2),fback);
   
    figure
    subplot(2,3,[1,4]);
    bode(H);
    title('system with feedback');
    subplot(2,3,[2,3]);
    pzmap(H);
    subplot(2,3,[5,6]);
    impulse(H);
    pole(H)
%% part 9(another system)

wn=100;
u=.3;
h2=tf(wn^2,[1,2*u*wn,wn^2]);

Tz=logspace(-5,100,20);
Tz=[0,Tz];
for i=Tz
    H=tf([i,1],1);
    figure
    subplot(2,3,[1,4]);
    bode(H);
    subplot(2,3,[2,3]);
    pzmap(H);
    subplot(2,3,[5,6]);
    impulse(H);
end

%%
Tz=2;
h1=tf([Tz,1],1);
H=series(h1,h2);
    figure
    subplot(2,3,[1,4]);
    bode(H);
    subplot(2,3,[2,3]);
    pzmap(H);
    subplot(2,3,[5,6]);
    impulse(H);

%% Question 3 part 1
T=.01;
a=[.2,.5,.8,1,1.5,2,3];

for i=a
    H=tf(T,[1,i-1],T);
    figure
    impulse(H)
    txt=sprintf(' = %f',i);
    title(['impulse response given \alpha ',txt])
end
%% part 2
T=.01;
a=[.2,1,2,3];

for i=a
    H=tf([T,0],[1+i,-1],T);
    figure
    impulse(H)
    txt=sprintf(' = %f',i);
    title(['impulse response given \alpha ',txt])
end
%% part 3
for i=a
    H=tf([T,0,0],[1.5+i,-2,.5],T);
    figure
    impulse(H)
    txt=sprintf(' = %f',i);
    title(['impulse response given \alpha ',txt])
end
%% 3 b
for i=a
    H=tf([T,0,0],[-1/12,2/3,i,-2/3,1/12],T);
    figure
    impulse(H)
    txt=sprintf(' = %f',i);
    title(['impulse response given \alpha ',txt])
end

