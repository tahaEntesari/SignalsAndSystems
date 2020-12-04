load("E:\books\Signal & System\EEGDATA\HW3\featu.mat");
%% calculate j-value for  persons ,var
for per=1:7 %persons
clear A B
Mu0_var = cal_Mu0(feature(2*per-1).var);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
              A.j_var.feature(i).bet(k,j) = cal_j(feature(2*per-1).var(k).Train_Data(i,:,:), feature(2*per-1).var(j).Train_Data(i,:,:), Mu0_var(i));
            end
        end 
     end
end
clear Mu0_var
%% calculate j-value for 3rd person ,mean
Mu0_mean = cal_Mu0(feature(2*per-1).mean);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_mean.feature(i).bet(k,j) = cal_j(feature(2*per-1).mean(k).Train_Data(i,:,:), feature(2*per-1).mean(j).Train_Data(i,:,:), Mu0_mean(i));
            end
        end 
     end
end
clear Mu0_mean
%% calculate j-value for 3rd person ,delta_band_energy

Mu0_delta_band_energy = cal_Mu0(feature(2*per-1).delta_band_energy);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_delta_band_energy.feature(i).bet(k,j) = cal_j(feature(2*per-1).delta_band_energy(k).Train_Data(i,:,:), feature(2*per-1).delta_band_energy(j).Train_Data(i,:,:), Mu0_delta_band_energy(i));
            end
        end 
     end
end
clear Mu0_delta_band_energy 
%% calculate j-value for 3rd person ,alpha_band_energy

Mu0_alpha_band_energy = cal_Mu0(feature(2*per-1).alpha_band_energy);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_alpha_band_energy.feature(i).bet(k,j) = cal_j(feature(2*per-1).alpha_band_energy(k).Train_Data(i,:,:), feature(2*per-1).alpha_band_energy(j).Train_Data(i,:,:), Mu0_alpha_band_energy(i));
            end
        end 
     end
end
clear Mu0_alpha_band_energy 

%% calculate j-value for 3rd person ,tetha_band_energy

Mu0_tetha_band_energy = cal_Mu0(feature(2*per-1).tetha_band_energy);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_tetha_band_energy.feature(i).bet(k,j) = cal_j(feature(2*per-1).tetha_band_energy(k).Train_Data(i,:,:), feature(2*per-1).tetha_band_energy(j).Train_Data(i,:,:), Mu0_tetha_band_energy(i));
            end
        end 
     end
end
clear Mu0_tetha_band_energy 
%% calculate j-value for 3rd person ,beta_band_energy

Mu0_beta_band_energy = cal_Mu0(feature(2*per-1).beta_band_energy);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_beta_band_energy.feature(i).bet(k,j) = cal_j(feature(2*per-1).beta_band_energy(k).Train_Data(i,:,:), feature(2*per-1).beta_band_energy(j).Train_Data(i,:,:), Mu0_beta_band_energy(i));
            end
        end 
     end
end
clear Mu0_beta_band_energy 
%% calculate j-value for 3rd person ,form_factor

Mu0_form_factor = cal_Mu0(feature(2*per-1).form_factor);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_form_factor.feature(i).bet(k,j) = cal_j(feature(2*per-1).form_factor(k).Train_Data(i,:,:), feature(2*per-1).form_factor(j).Train_Data(i,:,:), Mu0_form_factor(i));
            end
        end 
     end
end
clear Mu0_form_factor 
%% calculate j-value for 3rd person ,meanfreq

Mu0_meanfreq = cal_Mu0(feature(2*per-1).meanfreq);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_meanfreq.feature(i).bet(k,j) = cal_j(feature(2*per-1).meanfreq(k).Train_Data(i,:,:), feature(2*per-1).meanfreq(j).Train_Data(i,:,:), Mu0_meanfreq(i));
            end
        end 
     end
end
clear Mu0_meanfreq 
%% calculate j-value for 3rd person ,medfreq

Mu0_medfreq = cal_Mu0(feature(2*per-1).medfreq);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_medfreq.feature(i).bet(k,j) = cal_j(feature(2*per-1).medfreq(k).Train_Data(i,:,:), feature(2*per-1).medfreq(j).Train_Data(i,:,:), Mu0_medfreq(i));
            end
        end 
     end
end
clear Mu0_medfreq 
%% calculate j-value for 3rd person ,modfreq

Mu0_modefreq = cal_Mu0(feature(2*per-1).modefreq);
for k = 1:4
     for j=k+1:4
        for i = 1:64
            if (k~=j)
               A.j_modefreq.feature(i).bet(k,j) = cal_j(feature(2*per-1).modefreq(k).Train_Data(i,:,:), feature(2*per-1).modefreq(j).Train_Data(i,:,:), Mu0_modefreq(i));
            end
        end 
     end
end
clear Mu0_modefreq 

%% calculate j-value for 3rd person ,dwt
% calculate mu0
Mu0_dwt = zeros(64,183);
for i=1:64
    for j=1:183
        Mu0_dwt(i,j) =  (mean(feature(2*per-1).dwt(1).Train_Data(i,j)) + mean(feature(2*per-1).dwt(2).Train_Data(i,j)) + mean(feature(2*per-1).dwt(3).Train_Data(i,j)) + mean(feature(2*per-1).dwt(4).Train_Data(i,j)))/4;
    end
end
% calculte j-val
for m=1:4
    for n= m+1 : 4
        if (m~=n)
            for i = 1:64
                   for j=1:183
                        A.j_dwt.feature(i,j).bet(m,n) = cal_j(feature(2*per-1).dwt(m).Train_Data(i,j,:), feature(2*per-1).dwt(n).Train_Data(i,j,:), Mu0_dwt(i,j));
                   end 
            end
        end
    end
end
clear Mu0_dwt
%% calculate j-value for 3rd person , DST
% calculate mu0
Mu0_DST = zeros(64,360);
for i=1:64
    for j=1:183
        Mu0_DST(i,j) =  (mean(feature(2*per-1).DST(1).Train_Data(i,j)) + mean(feature(2*per-1).DST(2).Train_Data(i,j)) + mean(feature(2*per-1).DST(3).Train_Data(i,j)) + mean(feature(2*per-1).DST(4).Train_Data(i,j)))/4;
    end
end
% calculte j-val
for m=1:4
    for n= m+1 : 4
        if (m~=n)
            for i = 1:64
                   for j=1:360
                        A.j_DST.feature(i,j).bet(m,n) = cal_j(feature(2*per-1).DST(m).Train_Data(i,j,:), feature(2*per-1).DST(n).Train_Data(i,j,:), Mu0_DST(i,j));
                   end 
            end
        end
    end
end
clear Mu0_DST

%% calculate j-value for 3rd person , DCT
% calculate mu0
Mu0_DCT = zeros(64,360);
for i=1:64
    for j=1:360
        Mu0_DCT(i,j) =  (mean(feature(2*per-1).DCT(1).Train_Data(i,j)) + mean(feature(2*per-1).DCT(2).Train_Data(i,j)) + mean(feature(2*per-1).DCT(3).Train_Data(i,j)) + mean(feature(2*per-1).DCT(4).Train_Data(i,j)))/4;
    end
end
% calculte j-val
for m=1:4
    for n= m+1 : 4
        if (m~=n)
            for i = 1:64
                   for j=1:360
                        A.j_DCT.feature(i,j).bet(m,n) = cal_j(feature(2*per-1).DCT(m).Train_Data(i,j,:), feature(2*per-1).DCT(n).Train_Data(i,j,:), Mu0_DCT(i,j));
                   end 
            end
        end
    end
end
clear Mu0_DCT
%% get mean of j-valuse in A
for i=1:64
    B.j_var(i).feature = mean_2dmat(A.j_var.feature(i).bet);
    B.j_mean(i).feature = mean_2dmat(A.j_mean.feature(i).bet);
    B.j_delta_band_energy(i).feature = mean_2dmat(A.j_delta_band_energy.feature(i).bet);
    B.j_alpha_band_energy(i).feature = mean_2dmat(A.j_alpha_band_energy.feature(i).bet);
    B.j_tetha_band_energy(i).feature = mean_2dmat(A.j_tetha_band_energy.feature(i).bet);
    B.j_beta_band_energy(i).feature = mean_2dmat(A.j_beta_band_energy.feature(i).bet);
    B.j_form_factor(i).feature = mean_2dmat(A.j_form_factor.feature(i).bet);
    B.j_meanfreq(i).feature = mean_2dmat(A.j_meanfreq.feature(i).bet);
    B.j_medfreq(i).feature = mean_2dmat(A.j_medfreq.feature(i).bet);
    B.j_modefreq(i).feature = mean_2dmat(A.j_modefreq.feature(i).bet);
end
for i=1:64
    for j=1:183
        B.j_dwt.feature(i,j) = mean_2dmat(A.j_dwt.feature(i,j).bet);
    end
end
for i=1:64
    for j=1:360
        B.j_DST.feature(i,j) = mean_2dmat(A.j_DST.feature(i,j).bet);
    end
end
for i=1:64
    for j=1:360
        B.j_DCT.feature(i,j) = mean_2dmat(A.j_DCT.feature(i,j).bet);
    end
end
%% mean of all j and their std
c = [struct2cell(B.j_var), struct2cell(B.j_mean), struct2cell(B.j_delta_band_energy), struct2cell(B.j_tetha_band_energy), struct2cell(B.j_alpha_band_energy), struct2cell(B.j_beta_band_energy), struct2cell(B.j_meanfreq), struct2cell(B.j_medfreq), struct2cell(B.j_medfreq), struct2cell(B.j_modefreq)];
d = [struct2cell(B.j_dwt), struct2cell(B.j_DCT), struct2cell(B.j_DST)];
a = [reshape(d{1},1,64*183), reshape(d{2},1,64*360), reshape(d{3},1,64*360)];
k = [a, [c{:}]];
mean_j_all = mean(k);
std_j_all = std(k);
clear c d a k % temporary variables
%% select feaures
for i = 1:64
    for j=1:4
        if(B.j_var(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).var(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_mean(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).mean(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_delta_band_energy(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).delta_band_energy(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_tetha_band_energy(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).tetha_band_energy(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_alpha_band_energy(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).alpha_band_energy(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_beta_band_energy(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).beta_band_energy(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_form_factor(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).form_factor(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_meanfreq(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).meanfreq(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_modefreq(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).modefreq(j).Train_Data(i,1,:) = 0;
        end
        if(B.j_medfreq(i).feature < mean_j_all + 3*std_j_all)
            feature(2*per-1).medfreq(j).Train_Data(i,1,:) = 0;
        end 
    end
end
for i=1:64
    for j=1:183
        for k = 1:4
            if(B.j_dwt.feature(i,j) < mean_j_all + 3*std_j_all)
                feature(2*per-1).dwt(k).Train_Data(i,j,:) = 0;
            end
        end 
    end
end
for i=1:64
    for j=1:360
        for k = 1:4
            if(B.j_DST.feature(i,j) < mean_j_all + 5*std_j_all)
                feature(2*per-1).DST(k).Train_Data(i,j,:) = 0;
            end
            if(B.j_DCT.feature(i,j) < mean_j_all + 3*std_j_all)
                feature(2*per-1).DCT(k).Train_Data(i,j,:) = 0;
            end
        end 
    end
end
end

