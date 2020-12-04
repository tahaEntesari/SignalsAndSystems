clc
clear
%%z function myconv
%% a   
n=-100:1:100; %sampling
if ( 0<=n && 5>=n )
    x=1;
else
    x=0;
end
h1 = exp(-n).*(unistep(n) - unistep(n-5)); %?  %function unistep
y1 = MyConv(x1,h1) %function MyConv
y11 = conv(x1,h1); %Convolution
subplot(3,1,1);
stem(n,x1); %?
title('x1[n]')
subplot(3,1,2);
stem(n,h1);
title('h1[n]')
subplot(3,1,3);
stem(n,y1);
title('y1[n]')
hold on 
subplot(3,1,1);
stem(n,x1);
title('x1[n]')
subplot(3,1,2);
stem(n,h1);
title('h1[n]')
subplot(3,1,3);
stem(n,y11);
title('y11[n]')
%% b
t = linspace(-100,100,200000); %sampling
x2 = -1.*unistep(t+1) + 3.*unistep(t) - 2.*unistep(t-1); %function unistep
h2 = unistep(t) - unistep(t-1).*unistep(t) - unistep(t-1); %function unistep
y2 = MyConv(x2,h2) %function MyConv
y22 = conv(x2,h2); %Convolution
subplot(3,1,1);
plot(t,x2);
title('x2(t)')
subplot(3,1,2);
plot(t,h2);
title('h2(t)')
subplot(3,1,3);
plot(t,y2);
title('y2(t)')
hold on 
subplot(3,1,1);
plot(t,x2);
title('x2(t)')
subplot(3,1,2);
plot(t,h2);
title('h2(t)')
subplot(3,1,3);
plot(t,y22);
title('y22(t)')
%% c
t = linspace(-100,100,200000); %sampling
f1=1;
f2=100;
x3 = sin(2*pi*f1*t) + sin(2*pi*f2*t);
if ( t <= 10/f1 && t >= -10/f1 )
    h3 = sin(2*f1*t) ;
else 
    h3=0;
    end 
unistep(t) - unistep(t-1).*unistep(t) - unistep(t-1); %function unistep
y3 = MyConv(x3,h3) %function MyConv
y33 = conv(x3,h3); %Convolution
subplot(3,1,1);
plot(t,x3);
title('x3(t)')
subplot(3,1,2);
plot(t,h3);
title('h3(t)')
subplot(3,1,3);
plot(t,y3);
title('y3(t)')
hold on 
subplot(3,1,1);
plot(t,x3);
title('x3(t)')
subplot(3,1,2);
plot(t,h3);
title('h3(t)')
subplot(3,1,3);
plot(t,y33);
title('y33(t)')
%% d 
n=-100:1:100; %sampling
x4 = delta(n) + 2.*delta(n-4) - delta(n+4); %function delta
h4 = 0.5^(n).*(unistep(n) - unistep(n-1)); %function unistep
y4 = MyConv(x1,h1) %function MyConv
y44 = conv(x1,h1); %Convolution
subplot(3,1,1);
stem(n,x4); %?
title('x4[n]')
subplot(3,1,2);
stem(n,h4);
title('h4[n]')
subplot(3,1,3);
stem(n,y4);
title('y4[n]')
hold on 
subplot(3,1,1);
stem(n,x4);
title('x4[n]')
subplot(3,1,2);
stem(n,h4);
title('h4[n]')
subplot(3,1,3);
stem(n,y44);
title('y44[n]')
