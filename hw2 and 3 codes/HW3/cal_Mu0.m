function out = cal_Mu0(A)
out = zeros(1,64);

for i=1:64
      out(i) =  (mean(A(1).Train_Data(i,1)) + mean(A(2).Train_Data(i,1)) + mean(A(3).Train_Data(i,1)) + mean(A(4).Train_Data(i,1)))/4;
end

