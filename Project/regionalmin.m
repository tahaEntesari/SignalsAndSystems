function [out] = regionalmin(x,n)
%
out=zeros(1,length(x));
for i=1:length(x)
    if(i-n<1)
        start=1;
    else
        start=i-n;
    end
    if(i+n>length(x))
        endl=length(x);
    else
        endl=i+n;
    end
    if(all(x(i)<x([start:(i-1),(i+1):endl])))
        out(i)=1;
    end
end



end

