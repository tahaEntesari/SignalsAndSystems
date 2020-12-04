%% Q3
clc 
clear 
close all

im=imread('coins3.jpg');
im=rgb2gray(im);

im=im>40;
siz=size(im);
ed=edge(im);
edge=ed;

%{
im=im2double(im);


im(im<0.109)=0;
im(im>=0.109)=1;

%imshow(p3)

siz=size(im);
edge=zeros(siz(1),siz(2));
for i=2:(siz(1)-1)
    for j=2:(siz(2)-1)
       if (abs(im(i,j)-im(i,j-1))==1 && (im(i,j)==im(i,j+1)))||(abs(im(i,j)-im(i-1,j))==1 && (im(i,j)==im(i+1,j)))
           edge(i,j)=1;
       end
    end
end
%}

%%
cir_finder=zeros(siz(1),siz(2),200);
r=zeros(1,4);
for i=2:(siz(1)-1)
    for j=2:(siz(2)-1)
        if edge(i,j)==1
            r(1)=i-1;
            r(2)=j-1;
            r(3)=siz(1)-i;
            r(4)=siz(2)-j;
            for k=5:min(r)
                for l=0:360
                  a=round(j-k*sin(l*pi/180));
                  b=round(i-k*cos(l*pi/180));
                  cir_finder(b,a,k)=cir_finder(b,a,k)+1;  
                end
            end      
        end
    end
end

surf(cir_finder(:,:,20))


%%
[s1,s2,s3]=size(cir_finder);
cir=cir_finder;
for k=1:s3
    a=(mean(mean(cir_finder(:,:,k)))+1.5*mean(var(cir_finder(:,:,k))));
    for i=1:s1
        for j=1:s2
            if(cir_finder(i,j,k)<a)
                cir(i,j,k)=0;
            
            end
        end
    end
    
end

cir2=sum(cir(:,:,30:50),3);

cir2(cir2<100)=0;
%%
a=imregionalmax(cir2);
[b1,b2]=find(a==1);
[b1,I]=sort(b1,'ascend');
b2=b2(I);
b=[b1,b2];
circle=zeros(length(b1),2);
w=1;
circle(1,:)=[b1(1),b2(1)];
weight=cir2(b1(1),b2(1));
for i=2:length(b1)
    if ((abs(b(i,1)-circle(w,1))>20) || (abs(b(i,2)-circle(w,2))>20))
        weight=cir2(b1(i),b2(i));
        w=w+1;
        circle(w,:)=b(i,:);
    else
        circle(w,1)=(circle(w,1)*weight+b(i,1)*cir2(b1(i),b2(i)))/(weight+cir2(b1(i),b2(i)));
        circle(w,2)=(circle(w,2)*weight+b(i,2)*cir2(b1(i),b2(i)))/(weight+cir2(b1(i),b2(i)));
        weight=weight+cir2(b1(i),b2(i));
    end
end
%%
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
%%
groups=zeros(s1,s2,length(circle));
for i=1:s1
    for j=1:s2
        if(edge(i,j)==1)
            c=circle-[i,j];
            d=sqrt(c(:,1).^2+c(:,2).^2);
            [m,k]=min(d);
            groups(i,j,k)=m;
        end
    end
end
radius=zeros(1,length(circle));
for i=1:length(circle)
    [a,b]=find(groups(:,:,i)~=0);
    radius(i)=sum(sum(groups(a,b,i)))/nnz(groups(a,b,i));
end
ci=[circle(:,2),circle(:,1)];
imshow(im)
viscircles(ci,radius)





