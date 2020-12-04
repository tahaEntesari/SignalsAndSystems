function c=MyConv(u,v)
    %------------------
    syms x
    k=1:1:length(u);
    u1 = u(k).*x.^(k-1);
    a = sum(u1);
    %------------------
    n=1:1:length(v);
    v1 = v(n).*x.^(n-1);
    b = sum(v1);
    %------------------
    f = a.*b;
    c = coeffs(f);
end