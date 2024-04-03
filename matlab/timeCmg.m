load fromJulia.mat
la = load("fromJulia_la.mat").la;
b = load("fromJulia_b.mat").b;

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime

tic();
pfun = cmg_sdd(la);
% pfun = cmg_precondition(la);
bt = toc();

tic();
[x,flag,relres,iter] = pcg(la, b, tol, maxits, pfun);
st = toc();

save toJulia st bt iter relres startTime

exit
