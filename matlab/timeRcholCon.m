
% safe call to rchol
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

bt = 0;
st = 0;
iter = 0;
tic();

n = length(la);
x = zeros(n,1);

surplus = sum(la, 2);

% Get the connected components of la
a = -la;
a(1:size(a, 1) + 1:end) = 0;

G = graph(a);
cc = conncomp(G);

% Number of connected components
numComponents = max(cc);

% Display the connected components
for i = 1:numComponents
    componentIndices = find(cc == i);
    % fprintf('Connected Component %d: %s\n', i, mat2str(componentIndices));
    if max(surplus(componentIndices)) > 100 * eps
        % SDDM 
        % nCon = length(componentIndices);
        laCon = la(componentIndices, componentIndices);
        p = amd(laCon);
        laConperm = laCon(p,p);
        LCon = rchol(laConperm);

        bt = bt + toc();

        tic();
        [xs, flag, relres, iterCon] = pcg(laConperm, b(componentIndices(p)), tol /numComponents, maxits, LCon, LCon');
        x(componentIndices(p)) =xs;
        iter = iter + iterCon;
        st = st + toc();
    else
        % Lap
        nCon = length(componentIndices);
        laConSub = la(componentIndices(1:nCon-1), componentIndices(1:nCon-1));
        p = amd(laConSub);
        laConSubperm = laConSub(p,p);
        LConSub = rchol(laConSubperm);

        bt = bt + toc();
        [xs, flag, relres, iterCon] = pcg(laConSubperm, b(componentIndices(p)), tol / numComponents, maxits, LConSub, LConSub');
        x(componentIndices(p)) = xs;
        x(componentIndices) = x(componentIndices) - mean(x(componentIndices));
        iter = iter + iterCon;
        st = st + toc();
    end
end

relres = norm(la*x-b)/norm(b)

save toJulia st bt iter relres startTime

exit
