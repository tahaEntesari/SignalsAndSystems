function [class] = multisvm(train,class1,class2,class3,class4,x)
% the classes must be 1 for the good property and 0 else 
warning off
c=zeros(1,4);
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
    
warning on
end

