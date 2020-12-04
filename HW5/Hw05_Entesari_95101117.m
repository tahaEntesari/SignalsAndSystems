%% Question 1
clear
image=imread('17.jpg');
image=im2double(image);
s=size(image);
gaussiannoise=imnoise(image,'gaussian');
poissonnoise=image.*poissrnd(1,[s(1),s(2)]);
saltpeppernoise=imnoise(image,'salt & pepper');
specklenoise=imnoise(image,'speckle');
uniformnoise=image+0.35*rand(s(1),s(2),s(3));
figure
subplot(3,2,1)
imshow(image);
title('Original image');
subplot(3,2,2)
imshow(uniformnoise);
title('Uniform additive noise');
subplot(3,2,3)
imshow(gaussiannoise);
title('Gaussian noise');
subplot(3,2,4)
imshow(poissonnoise);
title('Poisson multiplicative noise');
subplot(3,2,5)
imshow(saltpeppernoise);
title('Salt & pepper noise');
subplot(3,2,6)
imshow(specklenoise);
title('Speckle noise');
% filtering
% median filter

uniformmed=median2dfilter(uniformnoise,1,7);
gaussianmed=median2dfilter(gaussiannoise,1,7);
specklemed=median2dfilter(specklenoise,1,7);
saltpeppermed=median2dfilter(saltpeppernoise,1,7);
poissonmed=median2dfilter(poissonnoise,1,7);

% gaussian filter
uniformgauss=imgaussfilt(uniformnoise,3);
gaussiangauss=imgaussfilt(gaussiannoise,1.5);
specklegauss=imgaussfilt(specklenoise,2);
saltpeppergauss=imgaussfilt(saltpeppernoise,2);
poissongauss=imgaussfilt(poissonnoise,2);
% wiener filter
for i=1:3
    uniformwiener(:,:,i)=wiener2(uniformnoise(:,:,i),[7,7]);
    gaussianwiener(:,:,i)=wiener2(gaussiannoise(:,:,i),[10,10]);
    specklewiener(:,:,i)=wiener2(specklenoise(:,:,i),[20,20]);
    saltpepperwiener(:,:,i)=wiener2(saltpeppernoise(:,:,i),[25,25]);
    poissonwiener(:,:,i)=wiener2(poissonnoise(:,:,i),[20,20]);
end

% plotting
figure
subplot(2,2,1)
imshow(gaussiannoise);
title('Gaussian noise')
subplot(2,2,2)
imshow(gaussianmed);
title('Median filter applied')
subplot(2,2,3)
imshow(gaussiangauss);
title('Gaussian filter applied')
subplot(2,2,4)
imshow(gaussianwiener);
title('Wiener filter applied')

figure
subplot(2,2,1)
imshow(poissonnoise);
title('Poisson noise')
subplot(2,2,2)
imshow(poissonmed);
title('Median filter applied')
subplot(2,2,3)
imshow(poissongauss);
title('Gaussian filter applied')
subplot(2,2,4)
imshow(poissonwiener);
title('Wiener filter applied')

figure
subplot(2,2,1)
imshow(saltpeppernoise);
title('Salt & pepper noise')
subplot(2,2,2)
imshow(saltpeppermed);
title('Median filter applied')
subplot(2,2,3)
imshow(saltpeppergauss);
title('Gaussian filter applied')
subplot(2,2,4)
imshow(saltpepperwiener);
title('Wiener filter applied')

figure
subplot(2,2,1)
imshow(specklenoise);
title('Speckle noise')
subplot(2,2,2)
imshow(specklemed);
title('Median filter applied')
subplot(2,2,3)
imshow(specklegauss);
title('Gaussian filter applied')
subplot(2,2,4)
imshow(specklewiener);
title('Wiener filter applied')

figure
subplot(2,2,1)
imshow(uniformnoise);
title('Uniform additive noise')
subplot(2,2,2)
imshow(uniformmed);
title('Median filter applied')
subplot(2,2,3)
imshow(uniformgauss);
title('Gaussian filter applied')
subplot(2,2,4)
imshow(uniformwiener);
title('Wiener filter applied')


%% Question 2
clear
[y_sing,Fs_sing]=audioread('singing.wav');
[y_lord,Fs_lord]=audioread('sound.wav');
% interpolating data to have the same sampling freq
% upsampling
Fs_max_lcm=lcm(Fs_lord,Fs_sing);
t_lord_reformed=0:(1/Fs_max_lcm):(length(y_lord)/Fs_lord);
t_sing_reformed=0:(1/Fs_max_lcm):(length(y_sing)/Fs_sing);
t_lord_orig=0:(1/Fs_lord):(length(y_lord)/Fs_lord);
t_sing_orig=0:(1/Fs_sing):(length(y_sing)/Fs_sing);
y_lord_up=spline(t_lord_orig(1:end-1)',y_lord,t_lord_reformed);
y_sing_up=spline(t_sing_orig(1:end-1)',y_sing,t_sing_reformed);
% downsampling
y_lord_down=y_lord_up(1:Fs_max_lcm/Fs_sing:end);
y_sing_down=y_sing_up(1:Fs_max_lcm/Fs_sing:end);
% converting the arrays to the same size
y_lord=y_lord_down(1:length(y_sing));
y_sing=y_sing_down(1:length(y_sing));
Fs=Fs_sing;
clear Fs_lord Fs_sing Fs_max_lcm t_lord_orig t_lord_reformed t_sing_orig t_sing_reformed
clear y_lord_down y_lord_up y_sing_down y_sing_up
%% main signal
sound(y_sing,Fs)
%%  phase deleted
sound(ifft(abs(fft(y_sing))),Fs);
%% real part only
sound(ifft(real(fft(y_sing))),Fs);
%%  main signal
sound(y_lord,Fs)
%%  phase deleted
sound(ifft(abs(fft(y_lord))),Fs);
%% phase shift by pi/4
y=fft(y_lord);
sound(real(ifft(abs(y).*exp(1j*(angle(y)+pi/4)))),Fs);
%% phase reversal
y=fft(y_lord);
sound(real(ifft(abs(y).*exp(-1j*(angle(y))))),Fs);
%% real part only
sound(ifft(real(fft(y_lord))),Fs);
%% random amplitude multiplication
ran=rand(1,length(y_lord));
y=y_lord.*ran;
sound(y,Fs)
%% generating sounds with their own amplitude but different phase
sound1=abs(fft(y_sing)).*exp(1j*angle((fft(y_lord))));
sound2=abs(fft(y_lord)).*exp(1j*angle(fft(y_sing)));
%%
sound(ifft(sound1),Fs);
%%
sound(ifft(sound2),Fs);
%% Images
clear
im1=imread('01.jpg');
im2=imread('02.jpg');
fft_im1=fft(im1);
fft_im2=fft(im2);
abs1phase2=abs(fft_im1).*exp(1j*angle(fft_im2));
abs2phase1=abs(fft_im2).*exp(1j*angle(fft_im1));
imwrite(1/255*real(ifft(abs1phase2)),'03.jpg');
imwrite(1/255*real(ifft(abs2phase1)),'04.jpg');
figure
subplot(2,2,1);
imshow(1/255*(real((ifft(abs1phase2)))));
title('Normalized image 01 with phase of image 02');
subplot(2,2,2);
imshow(1/255*(real(ifft(abs2phase1))));
title('Normalized image 02 with phase of image 01');
subplot(2,2,3);
imshow(real(ifft(abs1phase2)));
title('image 01 with phase of image 02');
subplot(2,2,4);
imshow(real(ifft(abs2phase1)));
title('image 02 with phase of image 01');
%% Question 3
clear
im=imread('coins3.jpg');
im=im2double(im);
fft2plot(im,0)
title('fft absolute plot');
%% coin finding using matlab's own functions
siz=zeros(1,300);
for i=10:150
im=imread('coins3.jpg');
im=rgb2gray(im);
im=im>i;
stats = regionprops('table',im,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
siz(i)=length(radii);
end
siz(siz==0)=[];
fprintf('There are %i coins in this photo\n',min(siz));
%%  coin finding using hough transform
% loading the data and finding the edges
clear
im=imread('coins3.jpg');
rgb=im; %don't delete
im=rgb2gray(im);
im=im>40;
im2=edge(im);
[s1,s2]=size(im);
% worst case coin maximum radius estimation
indx1=zeros(1,s1);
indx2=zeros(1,s1);
for i=1:s1
    a=find(im2(i,:)==1);
    if(isempty(a))
        m=0;
        M=0;
    else
        m=min(a);
        M=max(a);
    end
    indx1(i)=m;
    indx2(i)=M;
end
indy1=zeros(1,s1);
indy2=zeros(1,s1);
for i=1:s2
    a=find(im2(:,i)==1);
    if(isempty(a))
        m=0;
        M=0;
    else
        m=min(a);
        M=max(a);
    end
    indy1(i)=m;
    indy2(i)=M;
end
r1=indx2-indx1;
r2=indy2-indy1;
c=[r1,r2];
rmax=sum(sum(([r1,r2])))/(2*nnz([r1,r2]));
% finding possible circle centers for the given circumferences
s3=round(max([s1,s2])/2);
circle=zeros(s1,s2,s3);
r=zeros(1,4);
for i=2:(s1-1)
    for j=2:(s2-1)
        if im2(i,j)==1
            r(1)=i-1;
            r(2)=j-1;
            r(3)=s1-i;
            r(4)=s2-j;
            for k=5:min(r)
                for l=0:360
                  a=round(j-k*sin(l*pi/180));
                  b=round(i-k*cos(l*pi/180));
                  circle(b,a,k)=circle(b,a,k)+1;  
                end
            end      
        end
    end
end
figure
% setting bound for unwanted small data
for k=1:s3
    b=circle(:,:,k);
    b=b(:);
    a=(sum(sum(circle(:,:,k)))/(nnz(circle(:,:,k))+1)+10*std(b));
    for i=1:s1
        for j=1:s2
            if(circle(i,j,k)<a)
                circle(i,j,k)=0;      
            end
        end
    end
end
circle2=sum(circle(:,:,20:floor(rmax)),3);
% finding the location of the maximas and thus the center of the circles
a=imregionalmax(circle2);
[b1,b2]=find(a==1);
[b1,I]=sort(b1,'ascend');
b2=b2(I);
b=[b1,b2];
circle=zeros(length(b1),2);
w=1;
circle(1,:)=[b1(1),b2(1)];
weight=circle2(b1(1),b2(1));
% assuming that the circle are at least 20 pixels away from each other
for i=2:length(b1)
    if ((abs(b(i,1)-circle(w,1))>20) || (abs(b(i,2)-circle(w,2))>20))
        weight=circle2(b1(i),b2(i));
        w=w+1;
        circle(w,:)=b(i,:);
    else
        circle(w,1)=(circle(w,1)*weight+b(i,1)*circle2(b1(i),b2(i)))/(weight+circle2(b1(i),b2(i)));
        circle(w,2)=(circle(w,2)*weight+b(i,2)*circle2(b1(i),b2(i)))/(weight+circle2(b1(i),b2(i)));
        weight=weight+circle2(b1(i),b2(i));
    end
end
b1=circle(:,1);
b2=circle(:,2);
[b2,I]=sort(b2,'ascend');
b1=b1(I);
b=[b1,b2];
circle=zeros(length(b1),2);
w=1;
circle(1,:)=[b1(1),b2(1)];
for i=2:length(b1)
    if ((abs(b(i,1)-circle(w,1))>20) || (abs(b(i,2)-circle(w,2))>20))
    w=w+1;
    circle(w,:)=b(i,:);
    else
        circle(w,1)=.5*(circle(w,1)+b(i,1));
        circle(w,2)=.5*(circle(w,2)+b(i,2));
    end
end
circle(circle(:,1)==0,:)=[];
% finding the radius of the circles
[s,~]=size(circle);
groups=zeros(s1,s2,s);
for i=1:s1
    for j=1:s2
        if(im2(i,j)==1)
            c=circle-[i,j];
            d=sqrt(c(:,1).^2+c(:,2).^2);
            [m,k]=min(d);
            groups(i,j,k)=m;
        end
    end
end
radius=zeros(s,1);
for i=1:s
    [a,b]=find(groups(:,:,i)~=0);
    radius(i)=sum(sum(groups(a,b,i)))/(nnz(groups(a,b,i)));
end
% merging circles that are close to each other

w=1;
ci2=[];

for i=1:s
    for j=(i+1):s
        distance=sqrt(sum((circle(i,:)-circle(j,:)).^2));
        if (distance<radius(i) ||distance<radius(j) )
            if(radius(i)>radius(j))
                ci2(w,:)=circle(j,:);
            else
                ci2(w,:)=circle(i,:);   
            end
            w=w+1;
        end
    end
end
[s,~]=size(ci2);
for i=1:s
	radius(find(circle==ci2(i,1)))=[];
    circle(find(circle==ci2(i,1)),:)=[];
end
% plotting the image and the circles in one figure
ci=[circle(:,2),circle(:,1)];
figure
imshow(rgb)
viscircles(ci,radius);
circles=[circle,radius];
fprintf([' the centers of the circles and the corresponding radius are:\n'...
   , 'Center Coordinates\tRadius\n']);
disp(circles)



