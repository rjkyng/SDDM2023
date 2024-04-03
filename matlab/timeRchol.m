
load fromJulia.mat

la = load("fromJulia_la.mat").la;
b = load("fromJulia_b.mat").b;

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime


% Add paths to rchol and rchol_lap
rcholHome = getenv('RCHOL_HOME');
addpath(fullfile(rcholHome));
addpath(fullfile(rcholHome, 'rchol'));
addpath(fullfile(rcholHome, 'rchol_lap'));

tic();
n = length(la)
p = amd(la);
laperm = la(p,p);
L = rchol(laperm);

bt = toc();

tic();
[xs,flag,relres,iter] = pcg(laperm, b(p), tol, maxits, L, L');
x = zeros(n,1);
x(p) = xs;
st = toc();

relres = norm(la*x-b)/norm(b)


save toJulia st bt iter relres startTime

exit
