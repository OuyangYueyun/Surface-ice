function[F]=OLR(T,n)
    %调用原函数
    sigma=5.6703e-8;
    if n==0
        F=sigma*T^4;
    %调用一阶导
    elseif n==1
        F=4*sigma*T^3;
    end
end