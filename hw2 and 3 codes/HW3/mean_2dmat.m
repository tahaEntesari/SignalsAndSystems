function out = mean_2dmat(A)
[m,n] = size(A);
out = 0;
for i=1:m
    for j=1:n
        out = out + (A(i,j) ~= 0)* A(i,j);
    end
end

out = out / 6;
