function [] = fft2plot(im,which)
% which==1 means that 'im' is in frequency domain
    if(which==1)
    pf1=fftshift(im);
    pf1=mat2gray(log(abs(pf1)+1));
    imshow(pf1);
    else
        im=fft2(im);
        pf1=fftshift(im);
        pf1=mat2gray(log(abs(pf1)+1));
        imshow(pf1);
    end
end

