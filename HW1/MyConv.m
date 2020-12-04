function [out] = MyConv(A,B)
%
C=convmtx(A,numel(B));
out=B*C;
end

