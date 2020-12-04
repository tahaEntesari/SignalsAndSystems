function [class] = multisvm(train,class1,class2,class3,class4,x,which)
% which =1 uses predict and fitcsvm
% which=0 uses svmclassify and svmtrain
% the classes must be 1 for the good property and 0 else 
warning off
c=zeros(1,4);
if which==1
    c(1)=predict(fitcsvm(train,class1),x);
    c(2)=predict(fitcsvm(train,class2),x);
    c(3)=predict(fitcsvm(train,class3),x);
    c(4)=predict(fitcsvm(train,class4),x);
    if(sum(c)>1) 
        class=NaN;
        fprintf('predict svm grouped this into more than 1 group\n');
    elseif sum(c)==0
        class=NaN;
        fprintf('predict svm did not group this into any group\n');
    else 
        class=(find(c==1));
        fprintf('predict svm classified this into group %d\n',class);
    end
    
elseif which==0
    c(1)=svmclassify((svmtrain(train,class1)),x);
    c(2)=svmclassify((svmtrain(train,class2)),x);
    c(3)=svmclassify((svmtrain(train,class3)),x);
    c(4)=svmclassify((svmtrain(train,class4)),x);
  
    if(sum(c)>1) 
        class=NaN;
  %      fprintf('svmclassify grouped this into more than 1 group\n');
    elseif sum(c)==0
        class=NaN;
     %   fprintf('svmclassify did not group this into any group\n');
    else 
        class=(find(c==1));
       % fprintf('svmclassify classified this into group %d\n',class);
    end
    %}
end
warning on
end

