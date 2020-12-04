% til about line 110 is only loading,filtering and downsampling.
%no need for these to be done again.
%one can begin the code from the part "load newdata.mat"

%% loading the data

time=3;
path='E:\books\Signal & System\EEGDATA\HW3\data';
b=dir(path);

for i=3:length(b)
    txt=[b(i).folder,'\',b(i).name];
    data(i-2).name=b(i).name;
    if(length(strfind(b(i).name,'Execution'))==0)
        data(i-2).isexe=0;
    else
        data(i-2).isexe=1;
    end
    load(txt)
    data(i-2).Test=Test_Data_Scrambled;
    data(i-2).Train=Train_Data;
    %data(i-2).Test=load(txt,'Test_Data_Scrambled');
    %data(i-2).Train=load(txt,'Train_Data');
end
clear b path txt i Test_Data_Scrambled Train_Data

%% reducing mean from signals

for i=1:length(data)
    data(i).Test=data(i).Test-mean(data(i).Test(:,:,:),2);
    for j=1:length(data(i).Train)
        data(i).Train{j}=data(i).Train{j}-mean(data(i).Train{j}(:,:,:),2);
  
    end
end
clear i j


%%  filtering the sign

for i=1:length(data)
    for k=1:64
        [~,~,t]=size(data(i).Test);
        for m=1:t
          data(i).Test(k,:,m)=lowfilter(data(i).Test(k,:,m));  
        end
    end
    for j=1:length(data(i).Train)
        for k=1:64
            [~,~,t]=size(data(i).Train{j});
            for m=1:t
              data(i).Train{j}(k,:,m)=lowfilter(data(i).Train{j}(k,:,m));  
            end
        end
    end
end
clear i j k m t
%% saving 
save('data.mat','data','-v7.3')

%% downsampling
Fs=120;
for i=1:length(data)
    newdata(i).name=data(i).name;
    newdata(i).isexe=data(i).isexe;
    for k=1:64
        [~,~,t]=size(data(i).Test);
        
        for m=1:t
            newdata(i).Test(k,:,m)=data(i).Test(k,1:20:end,m);
            
        end
    end
    for j=1:length(data(i).Train)
        for k=1:64
            [~,~,t]=size(data(i).Train{j});
            for m=1:t
              newdata(i).Train{j}(k,:,m)=data(i).Train{j}(k,1:20:end,m);
              
            end
        end
    end
end
clear i j k m t data
%%
save('newdata.mat','newdata','Fs','time');
%%
load newdata.mat
%% low filter. pass band 45
for i=1:length(newdata)
    for k=1:64
        [~,~,t]=size(newdata(i).Test);
        for m=1:t
          newdata(i).Test(k,:,m)=lowfilter(newdata(i).Test(k,:,m));  
        end
    end
    for j=1:length(newdata(i).Train)
        for k=1:64
            [~,~,t]=size(newdata(i).Train{j});
            for m=1:t
              newdata(i).Train{j}(k,:,m)=lowfilter(newdata(i).Train{j}(k,:,m));  
            end
        end
    end
end
clear i j k m t
%% highpass filter passband .5hertz
for i=1:length(newdata)
    for k=1:64
        [~,~,t]=size(newdata(i).Test);
        for m=1:t
          newdata(i).Test(k,:,m)=highfilter(newdata(i).Test(k,:,m));  
        end
    end
    for j=1:length(newdata(i).Train)
        for k=1:64
            [~,~,t]=size(newdata(i).Train{j});
            for m=1:t
              newdata(i).Train{j}(k,:,m)=highfilter(newdata(i).Train{j}(k,:,m));  
            end
        end
    end
end
clear i j k m t
%% 
for i=1:14
    index=strfind(newdata(i).name,'_');
    newdata(i).subjectnumber=str2num(newdata(i).name((index(1)+1):(index(2)-1)));
end
clear i index
%%
L=time*Fs;
f = Fs*(0:(L/2))/L;
%% sth in the middle
save('newdata3.mat','newdata')
%% highpassed lowdoublepasses signal save
save('newdata2.mat','newdata','Fs','time','f','L');
%%  
load newdata2.mat
%%  calculating var and skewness
k=1;
for i=1:2:length(newdata)
          feature(k).subjectnumber=newdata(i).subjectnumber;
    for j=1:length(newdata(i).Train)
          feature(k).variance(j).Train=var(newdata(i).Train{j}(:,:,:)); 
          feature(k).skewness(j).Train=skewness(newdata(i).Train{j}(:,:,:),1,2); 
    end
    k=k+1;
end
clear i j k
%% calculating correlation between all channels for an activity one at a time for Trains
%calculating correlation between all activities for a channel one at a time for Trains
w=1;
for i=1:2:length(newdata)
    for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:t
            feature(w).channelcorrelation(j).testnumber(k).act=corr(newdata(i).Train{j}(:,:,k)');
        end
        for k=1:64
                feature(w).testcorrelation(j).channelnumber(k).act=...
                    corr(permute(newdata(i).Train{j}(k,:,:),[2,3,1]));
        end
    end
    w=w+1;
end
clear i j k t w
%% draw histogram
close all
k=13;
l=4;
channel=1;
[~,~,t]=size(newdata(k).Train{l});
for i=1:t
    figure;
    histogram(abs(newdata(k).Train{l}(channel,:,i)))
    txt=sprintf('subject %d, channel %d, class %d, test number %d'...
        ,newdata(k).subjectnumber,channel,l,i);
    title(txt);
    xlabel('amplitude');
    ylabel('density');
end
clear k l channel t i txt
%% Form factor
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                feature(w).formfactor(j).Train(k,l)=...
                    rms(newdata(i).Train{j}(k,:,l))/mean(abs(newdata(i).Train{j}(k,:,l)));
            end
        end
     end
    w=w+1;
end
clear i j k l t w
%% modefreq meanfreq medfreq
w=1;
for i=1:2:length(newdata)
    for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                P2 = abs(fft(newdata(i).Train{j}(k,:,l))/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                feature(w).modefreq(j).Train(k,l)=f(find(P1==max(P1)));  
                feature(w).meanfreq(j).Train(k,l)=meanfreq(newdata(i).Train{j}(k,:,l),Fs); 
                feature(w).medfreq(j).Train(k,l)=medfreq(newdata(i).Train{j}(k,:,l),Fs); 
            
            end
        end
    end
    w=w+1;
end
clear i j k w t l P1 P2
%%  dst dct
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                feature(w).dst(j).Train(k).act(l).channel=dst(newdata(i).Train{j}(k,:,l));
                feature(w).dct(j).Train(k).act(l).channel=dct(newdata(i).Train{j}(k,:,l));
            end
        end
     end
    w=w+1;
end
clear i j k l t w
%% band energy
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                feature(w).deltaenergy(j).Train(l,k)=...
                    fftenergy(deltafilter(newdata(i).Train{j}(k,:,l)),0,1);
            end
        end
     end
    w=w+1;
end
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
               feature(w).thetaenergy(j).Train(l,k)=...
                   fftenergy(thetafilter(newdata(i).Train{j}(k,:,l)),0,1);
            end
        end
     end
    w=w+1;
end
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
               feature(w).alphaenergy(j).Train(l,k)=...
                   fftenergy(alphafilter(newdata(i).Train{j}(k,:,l)),0,1);
            end
        end
     end
    w=w+1;
end
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                feature(w).betaenergy(j).Train(l,k)=...
                    fftenergy(betafilter(newdata(i).Train{j}(k,:,l)),0,1);
            end
        end
     end
    w=w+1;
end


clear i j k l t w
%%  dwt
w=1;
for i=1:2:14
     for j=1:length(newdata(i).Train)
        [~,~,t]=size(newdata(i).Train{j});
        for k=1:64
            for l=1:t
                feature(w).Dwt(j).Train(k,:,l)=dwt(newdata(w).Train{j}(k,:,l),'db1');
            end
        end
     end
    w=w+1;
end
clear i j k l t w
%%
save('feature.mat','feature');
%%
load newdata2.mat
load feature.mat
clear f Fs L time

%% j value
tic
for i=1:7
% finding mean of each feature.
    u0var=[];
    u0skew=[];
    u0form=[];
    
    u0medf=[];
    u0modf=[];
    u0meanf=[];
    
    u0dst=[];
    u0dct=[];
    u0dwt=[];
    
    u0beta=[];
    u0alpha=[];
    u0delta=[];
    u0theta=[];
    
    for j=1:4
        a=struct2cell(feature(i).variance(j));
        u0var=[u0var,mean([mean(a{1}(:,:,1:20))])];
        
        a=struct2cell(feature(i).skewness(j));
        u0skew=[u0skew,mean([mean(a{1}(:,:,1:20))])];
     
        a=struct2cell(feature(i).formfactor(j));
        u0form=[u0form,mean([mean(a{1}(:,1:20))])];
        
        a=struct2cell(feature(i).modefreq(j));
        u0modf=[u0modf,mean([mean(a{1}(:,1:20))])];
        
        a=struct2cell(feature(i).meanfreq(j));
        u0meanf=[u0meanf,mean([mean(a{1}(:,1:20))])];
        
        a=struct2cell(feature(i).medfreq(j));
        u0medf=[u0medf,mean([mean(a{1}(:,1:20))])];
        
        u0dwt=[u0dwt,mean(feature(i).dwt(j).Train_Data(:,:,1:20),3)];
        
        for k=1:20
            for l=1:64
        a=struct2cell(feature(i).dst(j).Train(l).act(k));
        u0dst=[u0dst,mean(a{:})];
        
        a=struct2cell(feature(i).dct(j).Train(l).act(k));
        u0dct=[u0dct,mean(a{:})];
        
        a=(feature(i).deltaenergy(j).Train(:,:));
        a=[a(:)];
        u0delta=[u0delta,mean(a)];
        
        a=(feature(i).alphaenergy(j).Train(:,:));
        a=[a(:)];
        u0alpha=[u0alpha,mean(a)];
        
        a=(feature(i).betaenergy(j).Train(:,:));
        a=[a(:)];
        u0beta=[u0beta,mean(a)];
        
        a=(feature(i).thetaenergy(j).Train(:,:));
        a=[a(:)];
        u0theta=[u0theta,mean(a)];
        
            end
        end    
    end
    u0skew=mean(u0skew);
    u0var=mean(u0var);
    u0form=mean(u0form);
    u0medf=mean(u0medf);
    u0modf=mean(u0modf);
    u0meanf=mean(u0meanf);
    u0dct=mean(u0dct);
    u0dst=mean(u0dst);
    u0theta=mean(u0theta);
    u0beta=mean(u0beta);
    u0alpha=mean(u0alpha);
    u0delta=mean(u0delta);
    u0dwt=mean(mean(u0dwt));
    for k=1:64
                for m=1:6
                    switch m  
                        case 1
                            n=1;l=2;
                        case 2
                        n=1;l=3;
                        case 3
                        n=1;l=4;
                        case 4
                        n=2;l=3;
                        case 5
                        n=2;l=4;
                        otherwise
                            n=3;l=4;
                    end
                    jval(i).variance(k,m)=jvalue(feature(i).variance(n).Train(k,:,:),feature(i).variance(l).Train(k,:,:),u0var);
                    jval(i).skewness(k,m)=jvalue(feature(i).skewness(n).Train(k,:,:),feature(i).skewness(l).Train(k,:,:),u0skew);
                    jval(i).formfactor(k,m)=jvalue(feature(i).formfactor(n).Train(k,:),feature(i).formfactor(l).Train(k,:),u0form);
                    jval(i).modefreq(k,m)=jvalue(feature(i).modefreq(n).Train(k,:),feature(i).modefreq(l).Train(k,:),u0modf);
                    jval(i).meanfreq(k,m)=jvalue(feature(i).meanfreq(n).Train(k,:),feature(i).meanfreq(l).Train(k,:),u0meanf);
                    jval(i).medfreq(k,m)=jvalue(feature(i).medfreq(n).Train(k,:),feature(i).medfreq(l).Train(k,:),u0medf);
                    
                    jval(i).alphaenergy(k,m)=jvalue(feature(i).alphaenergy(n).Train(:,k),feature(i).alphaenergy(l).Train(:,k),u0alpha);
                    jval(i).betaenergy(k,m)=jvalue(feature(i).betaenergy(n).Train(:,k),feature(i).betaenergy(l).Train(:,k),u0alpha);
                    jval(i).deltaenergy(k,m)=jvalue(feature(i).deltaenergy(n).Train(:,k),feature(i).deltaenergy(l).Train(:,k),u0alpha);
                    jval(i).thetaenergy(k,m)=jvalue(feature(i).thetaenergy(n).Train(:,k),feature(i).thetaenergy(l).Train(:,k),u0alpha);
                    for o=1:182
                        a=[feature(i).dwt(n).Train_Data(k,o,1:20)];
                        b=[feature(i).dwt(l).Train_Data(k,o,1:20)];
                        jval(i).dwt(k,o,m)=jvalue(a(:),b(:),u0dwt);
                    end
                    %{ 
                    for o=1:360
                        a=[];b=[];
                        for y=1:20
                            a=[a,feature(i).dst(n).Train(k).act(y).channel(o)];
                            b=[b,feature(i).dst(l).Train(k).act(y).channel(o)];
                        end
                            jval(i).dst(k,o,m)=jvalue(a,b,u0dst);
                        a=[];b=[];
                        for y=1:20
                            a=[a,feature(i).dct(n).Train(k).act(y).channel(o)];
                            b=[b,feature(i).dct(l).Train(k).act(y).channel(o)];
                        end
                            jval(i).dct(k,o,m)=jvalue(a,b,u0dct);
                    end
                    %}
                end
        end
        

end
clear u0alpha u0beta u0dct u0delta u0dst u0form u0meanf u0medf...
    u0modf u0skew u0theta u0var u0dwt i j k l m n a o b y


%%
for i=1:7
    jval(i).variance=mean(jval(i).variance,2);
    jval(i).skewness=mean(jval(i).skewness,2);
    jval(i).formfactor=mean(jval(i).formfactor,2);
    jval(i).modefreq=mean(jval(i).modefreq,2);
    jval(i).meanfreq=mean(jval(i).meanfreq,2);
    jval(i).medfreq=mean(jval(i).medfreq,2);
    jval(i).alphaenergy=mean(jval(i).alphaenergy,2);
    jval(i).betaenergy=mean(jval(i).betaenergy,2);
    jval(i).deltaenergy=mean(jval(i).deltaenergy,2);
    jval(i).thetaenergy=mean(jval(i).thetaenergy,2);
    jval(i).dst=mean(jval(i).dst,3);
    jval(i).dct=mean(jval(i).dct,3);
    jval(i).dwt=mean(jval(i).dwt,3);
end

%%
save('jval.mat','jval')
%%
save('newfeat.mat','feature')
%%
for i=1:7
    jval(i).tot=[jval(i).variance',jval(i).skewness',jval(i).formfactor'...
        jval(i).modefreq',jval(i).meanfreq',jval(i).medfreq'...
        jval(i).deltaenergy',jval(i).thetaenergy',jval(i).alphaenergy'...
        jval(i).betaenergy',reshape(jval(i).dst,1,64*360)...
        reshape(jval(i).dct,1,64*360),reshape(jval(i).dwt,1,64*182)];
end
feature(find(jval<mean(jval)+3*std(jval)))=0;


%% svm matrix 
for i=1:14
    for j=1:4
        
    a=struct2cell(feature(i).variance(j));
    a=permute(a{:},[3,1,2]);
    svm.variance(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
    a=struct2cell(feature(i).skewness(j));
    a=permute(a{:},[3,1,2]);
    svm.skewness(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
    a=struct2cell(feature(i).formfactor(j));
    a=permute(a{:},[2,1]);
    svm.formfactor(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
    a=struct2cell(feature(i).modefreq(j));
    a=permute(a{:},[2,1]);
    svm.modefreq(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
    a=struct2cell(feature(i).meanfreq(j));
    a=permute(a{:},[2,1]);
    svm.meanfreq(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
    a=struct2cell(feature(i).medfreq(j));
    a=permute(a{:},[2,1]);
    svm.medfreq(i,(20*(j-1)+1):(20*j),:)=a(1:20,:);
    
        for k=1:20
            a=struct2cell(feature(i).deltaenergy(j).Train(k).act);
            svm.deltaenergy(i,(20*(j-1)+k),:)=[a{:}];

            a=struct2cell(feature(i).thetaenergy(j).Train(k).act);
            svm.thetaenergy(i,(20*(j-1)+k),:)=[a{:}];

            a=struct2cell(feature(i).alphaenergy(j).Train(k).act);
            svm.alphaenergy(i,(20*(j-1)+k),:)=[a{:}];

            a=struct2cell(feature(i).betaenergy(j).Train(k).act);
            svm.betaenergy(i,(20*(j-1)+k),:)=[a{:}];
            
            a=struct2cell(feature(i).dst(j).Train(k).act);
            svm.dst(i,(20*(j-1)+k),:)=[a{:}];

            a=struct2cell(feature(i).dct(j).Train(k).act);
            svm.dct(i,(20*(j-1)+k),:)=[a{:}];
            
            a=(feature(i).dwt(j).Train_Data(:,:,k));
            svm.dwt(i,(20*(j-1)+k),:)=a(:)';
            
            
        end
        
    end
end
        svmmatrix=[permute(svm.variance,[1,3,2]),permute(svm.skewness,[1,3,2])...
            permute(svm.formfactor,[1,3,2]),permute(svm.modefreq,[1,3,2]),...
            permute(svm.meanfreq,[1,3,2]),permute(svm.medfreq,[1,3,2]),...
            permute(svm.deltaenergy,[1,3,2]),permute(svm.thetaenergy,[1,3,2]),...
            permute(svm.alphaenergy,[1,3,2]),permute(svm.betaenergy,[1,3,2]),...
            permute(svm.dst,[1,3,2]),permute(svm.dct,[1,3,2]),...
            permute(svm.dwt,[1,3,2])];
clear i j a svm
svmmatrix=permute(svmmatrix,[1,3,2]);
for i=1:7
    a=svmmatrix(i,1:15,:);
    b=svmmatrix(i,21:35,:);
    c=svmmatrix(i,41:55,:);
    d=svmmatrix(i,61:75,:);
    svm2train(i,:,:)=[a,b,c,d];
end
clear a b c d i
%%
save('svmorig.mat','svm2train','svmmatrix','svm_test');
%% svm 

class1=zeros(60,1);
class1(1:15)=1;
class2=zeros(60,1);
class2(16:30)=1;
class3=zeros(60,1);
class3(31:45)=1;
class4=zeros(60,1);
class4(46:60)=1;
[~,s1,s2]=size(svm2train);

b=svmmatrix(1,80,:);
b=b(:)';
multisvm(reshape(svm2train(1,:,:),s1,s2),class1,class2,class3,class4,b,0);

%% zeroing corresponding values of a_test

svm_test=a_test;
svmmatrix=a;

%%
[~,~,s]=size(svmmatrix);
tic
clear a
for j=1:7
    for i=1:s
        if((all(svmmatrix(j,:,i)==0)))
            svm_test(j,:,i)=0;   
        end
    end
end
%% deleting zeros
[~,~,s]=size(svmmatrix);
clear b c
for j=1:7
ii=1;
    for i=1:s
        if(~(all(svmmatrix(j,:,i)==0)))
            b(j,:,ii)=svmmatrix(j,:,i);
            c(j,:,ii)=svm_test(j,:,i);
            ii=ii+1;
        end
    end
end
svmmatrix=b;
svm_test=c;
clear b c

%%  kfold cross val for k=5
for mm=1:1
c=zeros(7,4);
vard=zeros(7);
for person=1:7
    for o=1:4
    u0=0;
    w=0;
    randed=rand(1,80);
    [b,index]=sort(randed);
    g=zeros(5,16);
    for i=1:5
        g(i,:)=index((16*(i-1)+1):16*i);
    end
    tic
    parfor i=1:5
        test=svmmatrix(person,g(i,:),:);

        switch i
            case 1
                 train=svmmatrix(person,g(2:5,:),:);
                 class=g(2:5,:);
                 class=[class(:)];

            case 2 
                train=svmmatrix(person,g([1,3:5],:),:);
                class=g([1,3:5],:);
                class=[class(:)];
            case 3
                train=svmmatrix(person,g([1,2,4,5],:),:);
                class=g([1,2,4,5],:);
                class=[class(:)];
            case 4
                train=svmmatrix(person,g([1:3,5],:),:);
                class=g([1:3,5],:);
                class=[class(:)];
            case 5
                train=svmmatrix(person,g(1:4,:),:);
                class=g(1:4,:);
                class=[class(:)];
        end

        classb=zeros(1,64);

                for j=1:64
                     if class(j)<21
                         classb(j)=1;
                     elseif class(j)>20&&class(j)<41
                         classb(j)=2;
                     elseif class(j)>40&&class(j)<61
                         classb(j)=3;
                     else
                         classb(j)=4;
                     end
                end

                classtest=g(i,:);
                classtest=[classtest(:)];
                classtestb=zeros(1,16);

                for j=1:16
                     if classtest(j)<21
                         classtestb(j)=1;
                     elseif classtest(j)>20&&classtest(j)<41
                         classtestb(j)=2;
                     elseif classtest(j)>40&&classtest(j)<61
                         classtestb(j)=3;
                     else
                         classtestb(j)=4;
                     end
                end


                class1=zeros(1,64);
                class1(find(classb==1))=1;

                class2=zeros(1,64);
                class2(find(classb==2))=1;

                class3=zeros(1,64);
                class3(find(classb==3))=1;

                class4=zeros(1,64);
                class4(find(classb==4))=1;
                for j=1:16
                    temp=test(1,j,:);
                    temp=[temp(:)]';

                    [~,~,s]=size(svmmatrix);
                    a=multisvm(reshape(train,64,s),class1,class2,class3,class4,temp,0);
                    if(a==classtestb(j))
                        w=w+1;
                        u0=u0+1;

                    end
                end

    end
    c(person,o)=sum(u0)*100/80;
    end
end
fprintf('for each person, the probability of correct answer is:\n');
disp(mean(c,2)/100);
vard=var(c,1,2);
fprintf('for each person, the variance of different cross validation\n probabilities of success is:\n');
disp(vard/100);
fprintf('the total probability of success:%f , error is %f\nmean variance of success: %f\n',mean(mean(c,2))/100,1-mean(mean(c,2))/100,mean(vard)/100);

end
%% svm on test

class1=zeros(1,80);
class2=zeros(1,80);
class3=zeros(1,80);
class4=zeros(1,80);
class1(1:20)=1;
class2(21:40)=1;
class3(41:60)=1;
class4(61:80)=1;

test_class=zeros(7,49);
[~,s1,s2]=size(svmmatrix);
parfor i=1:7
     for j=1:49
         temp=svm_test(i,j,:);
         temp=[temp(:)]';
         test_class(i,j)=multisvm(reshape(svmmatrix(i,:,:),s1,s2),class1,class2,class3,class4...
             ,temp,0);
     end
end



