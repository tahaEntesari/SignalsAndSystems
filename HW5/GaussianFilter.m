function out = GaussianFilter ( I , N , sigma )
% N = 5; 
% sigma = 2; 
% Generate Gaussian mask
ind = -floor(N/2) : floor(N/2);
[X,Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));

h = h / sum(h(:));
% Convert filter into a column vector
h = h(:);
% Filter our image
if(size(I,3)~=1)
    I=rgb2gray(I);
end
I = im2double(I);
I_pad = padarray(I, [floor(N/2) floor(N/2)]);
C = im2col(I_pad, [N N], 'sliding');
C_filter = sum(bsxfun(@times, C, h), 1);
out = col2im(C_filter, [N N], size(I_pad), 'sliding');
end
