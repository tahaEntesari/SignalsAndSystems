function [jval] = jvalue (x1,x2,u0)
%

u1=mean(x1);
u2=mean(x2);
s1=var(x1);
s2=var(x2);
jval=(abs(u0-u1)^2+abs(u0-u2)^2)/(s1^2+s2^2);
end

