function x= pcg_btv(x,A,T,mu)
  temp =T(:);
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
        y = A'*(A*x) + mu*diff1T(diff1(x));
    end
end