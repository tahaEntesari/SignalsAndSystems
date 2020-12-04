%% Question 1
clear

mm=sprintf('%d.jpg',16);
image=imread(mm);
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

%uniformmed=median2dfilter(uniformnoise,1,7);
%gaussianmed=median2dfilter(gaussiannoise,1,7);
%specklemed=median2dfilter(specklenoise,1,7);
%saltpeppermed=median2dfilter(saltpeppernoise,1,7);
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
