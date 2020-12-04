function [out]=median2dfilter(x,rgbORgray,n)
% rgbORgray=1 means that the image is rgb and thus the x is 3 dimentional
% rgbORgray=0 means gray scale.
% using square window of numofelements=9
% n, which must be odd, is the window square length

    [s1,s2,s3]=size(x);
    out=zeros(s1,s2,s3);
    out=squeeze(out);
    a=zeros(n,n);
    if (rgbORgray==1)
        y=zeros(3*s1,3*s2,3);
        y1=repmat(x(:,:,1),3);
        y2=repmat(x(:,:,2),3);
        y3=repmat(x(:,:,3),3);
        y(:,:,1)=y1;
        y(:,:,2)=y2;
        y(:,:,3)=y2;
        for i=1:3
            for j=(s1+1):(2*s1)
                for k=(s2+1):(2*s2)
                    for l=-floor(n/2):floor(n/2)
                        for m=-floor(n/2):floor(n/2)
                        a(l+ceil(n/2),m+ceil(n/2))=y(j+l,k+m,i);
                        end
                    end
                    out(j-s1,k-s2,i)=median(a(:));
                end
            end
        end
    else
        y=zeros(3*s1,3*s2);
        y=repmat(x,3);
       for j=(s1+1):(2*s1)
                for k=(s2+1):(2*s2)
                    for l=-floor(n/2):floor(n/2)
                        for m=-floor(n/2):floor(n/2)
                        a(l+ceil(n/2),m+ceil(n/2))=y(j+l,m+k);
                        end
                    end
                    out(j-s1,k-s2)=median(a(:));
                end
            end
    end
end