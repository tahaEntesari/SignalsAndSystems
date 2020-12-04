function out = cal_j (A,B,Mu0)

Mu1 = mean(A(1:20));
Mu2 = mean(B(1:20));
sigma1 = var(A(1:20));
sigma2 = var(B(1:20));

out = ((Mu0 - Mu1)^2 + (Mu0 - Mu2)^2)/(sigma1 + sigma2);