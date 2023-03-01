load fromJulia.mat

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime

tic();
n = length(la)
p = symrcm(la);
laperm = la(p,p);
L = ichol(laperm);

bt = toc();

tic();
[xs,flag,relres,iter] = pcg(laperm, b(p), tol, maxits, L, L');
x = zeros(n,1);
x(p) = xs;
st = toc();

relres = norm(la*x-b)/norm(b)


save toJulia st bt iter relres startTime

exit
