function y = myRelErr(x,y)
y = norm(x-y,2)/norm(x,2);
end