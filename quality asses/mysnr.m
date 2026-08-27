function y = mysnr(x0,x1)
y = 20*log10(norm(x0)/norm(x0-x1));