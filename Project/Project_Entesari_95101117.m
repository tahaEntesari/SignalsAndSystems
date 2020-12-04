%% Question 1
% importing the data
clear
cap_data=importdata("E:\books\Signal & System\Project\data - part1.xlsx");
[s1,s2]=size(cap_data);
capacitor_data=zeros(s1,s2);
% converting the data from cell format to double arrays
for i=1:s1
    for j=1:s2
        capacitor_data(i,j)=str2double(cap_data{i,j});
    end
end
clear i j cap_data
% deleting the first row of each run
capacitor_data(1:40:end,:)=[];
cap2=capacitor_data(:,1:12);
cap=zeros(8*30,12*39);
[s1,~]=size(capacitor_data);
% rearranging the matrix to only have 240 trials
for i=1:39
    w=1;
    for k=i:39:s1
      cap(w,(12*(i-1)+1):(12*(i))) = cap2(k,:);
      w=w+1;
    end
end
% labeling the new matrix
label=zeros(1,240);
for i=1:7
    label((30*i+1):(30*(i+1)))=i;
end
feature=[cap,label'];
% extra feature extraction 
% differentiation of consequent rows
dif=diff(cap2,1);
% the 39'th row is not meaningful since it measures the difference
% between the current run and the next one
dif(39:39:end,:)=[];
diff_feature=zeros(240,12*38);
for i=1:38
    w=1;
    for k=i:38:9120
      diff_feature(w,(12*(i-1)+1):(12*(i))) = dif(k,:);
      w=w+1;
    end
end

feature2=[cap,diff_feature,label'];
clear i k w

%%  Question 2
clear
chair_data=importdata("E:\books\Signal & System\Project\chair position dataset.xlsx");
[s1,~]=size(chair_data);
runs=zeros(1,21);
w=1;
for i=1:(s1-1)
    if(chair_data(i,11)~=chair_data(i+1,11))
    runs(w)=chair_data(i,12);
    w=w+1;
    end
    if(i==(s1-1))
        runs(w)=chair_data(i,12);
    end
end
chair_data(1:40:s1,:)=[];
chair2=chair_data(:,1:10);
chair=zeros(sum(runs),10*39);
[s1,~]=size(chair_data);
for i=1:39
    w=1;
    for k=i:39:s1
      chair(w,(10*(i-1)+1):(10*(i))) = chair2(k,:);
      w=w+1;
    end
end
label=zeros(1,sum(runs));
run=cumsum(runs);
for i=1:20
    label((run(i)+1):(run(i+1)))=i;
end
feature=[chair,label'];
feature(end,:)=[];
feature1=chair_data(:,1:11);
clear i k w
% supposing mean usage of data
for i=1:1222
    feature2(i,:)=mean(chair2((39*(i-1)+1):(39*i),:),1);
end
label2=label';
label2(end)=[];
feature2=[feature2,label2];
%% Question 3
clear
% importing deep breathing data
deep1=xlsread('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep1.csv','A:D');
deep2=xlsread('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep2.csv','A:D');
deep3=xlsread('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep3.csv','A:D');
deep4=zeros(400,6);
deep5=zeros(400,6);
deep4(:,1:4)=xlsread('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep4.csv','A:D');
deep5(:,1:4)=xlsread('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep5.csv','A:D');
deep11=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep1.csv');
deep11=deep11.textdata;
deep22=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep2.csv');
deep22=deep22.textdata;
deep33=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep3.csv');
deep33=deep33.textdata;
deep44=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep4.csv');
deep44=deep44.textdata;
deep55=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\deep5.csv');
deep55=deep55.textdata;
% extracting exact time
for i=1:400
    deep3(i,5)=str2double(deep33{i,6}(15:16));
    deep3(i,6)=str2double(deep33{i,6}(18:end));
    deep4(i,5)=str2double(deep44{i,6}(15:16));
    deep4(i,6)=str2double(deep44{i,6}(18:end));
    deep5(i,5)=str2double(deep55{i,6}(15:16));
    deep5(i,6)=str2double(deep55{i,6}(18:end));
end
for i=1:200
    deep1(i,5)=str2double(deep11{i,6}(15:16));
    deep1(i,6)=str2double(deep11{i,6}(18:end));
    deep2(i,5)=str2double(deep22{i,6}(15:16));
    deep2(i,6)=str2double(deep22{i,6}(18:end));
end
clear i deep11 deep22 deep33deep44 deep55 
deeptiming44=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\timing4.csv');
deeptiming55=importdata('E:\books\Signal & System\Project\Phase 3\Deep Breath\timing5.csv');
deeptiming4=zeros(8,8);
deeptiming5=zeros(9,8);
for i=1:8
    a=find(deeptiming44{i}==' ');
    a=a+4;
    b=a+3;
    for j=1:4
        deeptiming4(i,2*j-1)=str2double(deeptiming44{i}(a(j):(a(j)+1)));
        deeptiming4(i,2*j)=str2double(deeptiming44{i}(b(j):(b(j)+8)));
    end
end
for i=1:9
    a=find(deeptiming55{i}==' ');
    a=a+4;
    b=a+3;
    for j=1:4
        deeptiming5(i,2*j-1)=str2double(deeptiming55{i}(a(j):(a(j)+1)));
        deeptiming5(i,2*j)=str2double(deeptiming55{i}(b(j):(b(j)+8)));
    end
end
clear deeptiming44 deeptiming55 i j a b
% importing regular breathing data
regular1=xlsread('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular1.csv','A:D');
regular2=xlsread('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular2.csv','A:D');
regular3=zeros(400,6);
regular4=zeros(400,6);
regular3(:,1:4)=xlsread('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular3.csv','A:D');
regular4(:,1:4)=xlsread('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular4.csv','A:D');
reg1=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular1.csv');
reg2=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular2.csv');
reg3=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular3.csv');
reg4=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\regular4.csv');
reg1=reg1.textdata;
reg2=reg2.textdata;
reg3=reg3.textdata;
reg4=reg4.textdata;
% extracting exact time
for i=1:400
    regular1(i,5)=str2double(reg1{i,6}(15:16));
    regular1(i,6)=str2double(reg1{i,6}(18:end));
    regular2(i,5)=str2double(reg2{i,6}(15:16));
    regular2(i,6)=str2double(reg2{i,6}(18:end));
    regular3(i,5)=str2double(reg3{i,6}(15:16));
    regular3(i,6)=str2double(reg3{i,6}(18:end));
    regular4(i,5)=str2double(reg4{i,6}(15:16));
    regular4(i,6)=str2double(reg4{i,6}(18:end));
end
clear i reg4 reg3 reg2 reg1
%
regular3(1,1:4)=regular3(2,1:4);
regular2(1,1:4)=regular2(2,1:4);
timing3=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\timing3.csv');
timing4=importdata('E:\books\Signal & System\Project\Phase 3\Regular Breath\timing4.csv');
timing3start=timing3.data;
timing4start=timing4.data;
timing3end=timing3.textdata;
timing4end=timing4.textdata;
regulartiming3=zeros(20,2);
regulartiming4=zeros(24,2);
%
for i=1:20
    regulartiming3(i,3)=timing3start(i,1);
    regulartiming3(i,4)=timing3start(i,2);
    regulartiming3(i,1)=str2double(timing3end{i,2});
    regulartiming3(i,2)=str2double(timing3end{i,3}(1:9));
end
%
for i=1:24
    regulartiming4(i,3)=timing4start(i,1);
    regulartiming4(i,4)=timing4start(i,2);
    regulartiming4(i,1)=str2double(timing4end{i,2});
    regulartiming4(i,2)=str2double(timing4end{i,3}(1:9));
end
clear timing3 timing4 i timing3end timing4end timing3start timing4start
%% saving the data 
save data.mat
%% computation time saving
clear
load data.mat
% In this section and the next I will try to learn and see how the error of
% my algorithm is. in the section after these two I will print the results
% because at first I had wrote my code only for deep4 data, for testing on
% other data I only change the values of deep4 and deeptimingtiming4 with
% the desired values rather than changing the whole code
deep4=deep4;
deeptiming4=deeptiming4;
l2=length(deeptiming4);
l=length(deep4);
instart=zeros(1,l2);
inend=zeros(1,l2);
outstart=zeros(1,l2);
outend=zeros(1,l2);
threshold=0.2;
% finding index of startpoint of events
for i = 1:l2
    instart(i)=max(find(abs(deep4(:,6)-deeptiming4(i,2))<threshold));
    inend(i)=max(find(abs(deep4(:,6)-deeptiming4(i,4))<threshold));
    outstart(i)=max(find(abs(deep4(:,6)-deeptiming4(i,6))<threshold));
    b=max(find(abs(deep4(:,6)-deeptiming4(i,8))<threshold));
    if(isempty(b))
        outend(i)=l;
    elseif(i==l2 && b<outend(i-1))
        outend(i)=l;
    else
        outend(i)=b;
    end
end
clear i
dif=diff(deep4);
% the idea for this aproach yielded from the figure below
% plotyy(1:400,deep4(:,1),1:399,diff(deep4(:,1)));
l=length(deep4);
% regional maximums of the diff are where the inhaling starts
regmax=zeros(l-1,3);
% regional minimums of the diff are where the exhaling starts
regmin=zeros(l-1,3);
%dif(:,1)=filloutliers(dif(:,1),'spline');
neighboorhood=15;
for i=1:3
    regmax(:,i)=regionalmax(dif(:,i),neighboorhood);
    regmin(:,i)=regionalmin(dif(:,i),neighboorhood);
end

inhalestart=find(regmax(:,1)==1);
exhalestart=find(regmin(:,2)==1);
% inhale end is marked as the place where the "dif" passes zero
% likewise for exhale
if (exhalestart(1)<inhalestart(1)) 
    exhalestart(1)=[];
end 
% inhale/exhale end extraction
inhaleend=zeros(1,length(inhalestart));
exhaleend=l*ones(1,length(exhalestart));
for i=1:length(inhalestart)
    w=0;
    for k=inhalestart(i):(l-3)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            inhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
for i=1:length(exhalestart)
    w=0;
    for k=exhalestart(i):(l-3)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            exhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
clear i k l

%% regular breathing learning
clear
load data.mat
%
deep4=regular4;
deeptiming4=regulartiming4;
l2=length(deeptiming4);
instart=zeros(1,l2)';
outstart=zeros(1,l2)';
threshold=0.2;
for i = 1:(l2)
    instart(i)=max(find(abs(deep4(:,6)-deeptiming4(i,2))<threshold));
    outstart(i)=max(find(abs(deep4(:,6)-deeptiming4(i,2))<threshold));
end
clear i
%
dif=diff(deep4);
% the idea for this aproach yielded from the figure below
%plotyy(1:400,deep4(:,1),1:399,diff(deep4(:,1)));
l=length(deep4);
% regional maximums of the diff are where the inhaling starts
regmax=zeros(l-1,3);
% regional minimums of the diff are where the exhaling starts
regmin=zeros(l-1,3);
%dif(:,1)=filloutliers(dif(:,1),'spline');
neighboorhood=14;
for i=1:3
    regmax(:,i)=regionalmax(dif(:,i),neighboorhood);
    regmin(:,i)=regionalmin(dif(:,i),neighboorhood);
end
inhalestart=find(regmax(:,1)==1);
exhalestart=find(regmin(:,2)==1);
% inhale end is marked as the place where the "dif" passes zero
% likewise for exhale
if (exhalestart(1)<inhalestart(1)) 
    exhalestart(1)=[];
end 
% inhale/exhale end extraction
inhaleend=zeros(1,length(inhalestart))';
exhaleend=l*ones(1,length(exhalestart))';
for i=1:length(inhalestart)
    w=0;
    for k=inhalestart(i):(l-3)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            inhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
for i=1:length(exhalestart)
    w=0;
    for k=exhalestart(i):(l-4)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            exhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
clear i k l




%% On tests
clear
load data.mat
deep4=regular4;
l=length(deep4);
dif=diff(deep4);
% regional maximums of the diff are where the inhaling starts
regmax=zeros(l-1,3);
% regional minimums of the diff are where the exhaling starts
regmin=zeros(l-1,3);
neighboorhood=14;
for i=1:3
    regmax(:,i)=regionalmax(dif(:,i),neighboorhood);
    regmin(:,i)=regionalmin(dif(:,i),neighboorhood);
end

inhalestart=find(regmax(:,2)==1);
exhalestart=find(regmin(:,2)==1);
% inhale end is marked as the place where the "dif" passes zero
% likewise for exhale
if (exhalestart(1)<inhalestart(1)) 
    exhalestart(1)=[];
end 
% inhale/exhale end extraction
inhaleend=l*ones(1,length(inhalestart))';
exhaleend=l*ones(1,length(exhalestart))';
for i=1:length(inhalestart)
    w=0;
    for k=inhalestart(i):(l-4)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            inhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
for i=1:length(exhalestart)
    w=0;
    for k=exhalestart(i):(l-3)
        if (dif(k+1,1)*dif(k+2,1)<0 && dif(k+2,1)*dif(k+3,1)>0)
            exhaleend(i)=k+2;
            w=1;
        end
        if(w==1)
            break
        end
    end
end
clear i k 
% printing results
% inhaling
fprintf('The found results are as follows\n   Inhale start\t\tInhale end\n');
for i=1:length(inhalestart)
   disp([deep4(inhalestart(i),5),deep4(inhalestart(i),6),deep4(inhaleend(i),5),deep4(inhaleend(i),6)]);
    
end
% exhaling
fprintf('   Exhale start\t\tExhale end\n');
for i=1:length(exhalestart)
    disp([deep4(exhalestart(i),5),deep4(exhalestart(i),6),deep4(exhaleend(i),5),deep4(exhaleend(i),6)]);
end
clear i w

n=1:l;
inhale=zeros(1,l);
exhale=zeros(1,l);
for i=1:length(inhalestart)
    inhale(inhalestart(i):inhaleend(i))=1;
end
for i=1:length(exhalestart)
    exhale(exhalestart(i):exhaleend(i))=1;
end
close all
yyaxis right
plot(deep4(:,1:3));
yyaxis left
plot(n,inhale);
axis([1,l,0,5]);
hold on
plot(n,exhale,'r');
legend('inhale','exhale','channel 1','channel 2','channel 3')
