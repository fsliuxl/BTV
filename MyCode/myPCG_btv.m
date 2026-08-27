function x= myPCG_btv(x,A,T,N,mu)
  temp =T(:);
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
        x = reshape(x,N);
        y = A'*(A*x) + mu*diffT2(diff2(x,N),N);
        y = y(:);
    end
end