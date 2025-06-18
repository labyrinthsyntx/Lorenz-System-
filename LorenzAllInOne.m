%%========================================================================
%%  LORENZ ALL-IN-ONE SCRIPT
%%  DEREK MARTINEZ, ENGR.
%%========================================================================
clear; clc; close all; warning off backtrace;

%% 0 | USER PROMPT (baseline run)
pTxt   = {'sigma  (default 16):','beta  (default 4/3):','rho  (default 45):',...
          'tspan [t0 tf]:','IC [x0 y0 z0]:'};
dflts  = {'16','4/3','45','[0 100]','[1 1 1]'};
answ   = inputdlg(pTxt,'Baseline Lorenz parameters',1,dflts);
if isempty(answ), answ = dflts; end

sigma0 = str2double(answ{1});
beta0  = str2num(answ{2});    %#ok<ST2NM>
rho0   = str2double(answ{3});
tspan0 = str2num(answ{4});    %#ok<ST2NM>
x0     = str2num(answ{5})';   %#ok<ST2NM>
if numel(tspan0)~=2 || numel(x0)~=3
    error('Check tspan / initial-condition format.');
end

%% 1 | RESEARCH SUITE SETTINGS
dom.sigma      = [8 24];
dom.beta       = [1  4];
dom.rho        = [20 120];
nSamples       = 200;
nEnsembleIC    = 100;
bootstrapRuns  = 500;
tspanSweep     = [0 50];
dtSkip         = 5;
rng(1);

%% 2 | PARPOOL
if isempty(gcp('nocreate')), parpool; end

%% 3 | BASELINE RUN + PLOTS
fprintf('\n--- Baseline integration ----------------------------------\n');
sol0 = ode45(@(t,x)lorenzRHS(t,x,[sigma0 beta0 rho0]), tspan0, x0, ...
             odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1));
t0 = sol0.x;  X0 = sol0.y';

% 3-D phase
figure('Name','Baseline 3-D','Position',[100 100 820 600],'Color','w');
scatter3(X0(:,1),X0(:,2),X0(:,3),12,t0,'filled');
colormap jet; colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; title('Baseline attractor');

% Animation
figure('Name','Baseline animation','Color','w','Position',[80 80 800 600]);
h = plot3(X0(1,1),X0(1,2),X0(1,3),'b','LineWidth',1.4); hold on; grid on;
view(45,25); axis tight; xlabel X; ylabel Y; zlabel Z;
skip = max(floor(numel(t0)/800),1);
for i = 2:skip:numel(t0)
    set(h,'XData',X0(1:i,1),'YData',X0(1:i,2),'ZData',X0(1:i,3));
    drawnow limitrate;
end

% Projections
figure('Name','Baseline projections','Color','w','Position',[950 80 900 800]);
subplot(3,1,1), plot(X0(:,1),X0(:,2),'b'), title('XY'), grid on;
subplot(3,1,2), plot(X0(:,1),X0(:,3),'r'), title('XZ'), grid on;
subplot(3,1,3), plot(X0(:,2),X0(:,3),'g'), title('YZ'), grid on;

% Time series
figure('Name','Baseline series','Color','w','Position',[100 780 1000 620]);
plot(t0,X0,'LineWidth',1.2), legend('X','Y','Z'), grid on, title('State vs time');

% Lyapunov exponent
L0 = lyapExponent(X0(:,1), t0);
if isnan(L0)
    warning('lyapExponent returned NaN; using two-trajectory fallback.');
    eps0 = 1e-6;
    x1  = x0 + eps0*[1;0;0];
    sol1 = ode45(@(t,x)lorenzRHS(t,x,[sigma0 beta0 rho0]), tspan0, x1, ...
                 odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1));
    t1  = sol1.x;  X1 = sol1.y';
    X1i = interp1(t1, X1, t0);
    d0   = norm(x1 - x0);
    dEnd = norm(X1i(end,:) - X0(end,:));
    L0   = (1/(t0(end)-t0(1))) * log(dEnd/d0);
end
fprintf('Baseline Lyap = %.4f', L0);

% Correlation dimension (optional)
try
    fd0 = corrDimension(X0(:,1));
    fprintf(', CorrDim ≈ %.2f\n', fd0);
catch
    warning('corrDimension not found; skipping.');
    fprintf('\n');
end

%% 4 | LHS PARAMETER SWEEP
fprintf('\n--- LHS sweep (%d pts) -----------------------------------\n',nSamples);
lhs = lhsdesign(nSamples,3,'criterion','maximin','iterations',1e3);
sigS = dom.sigma(1) + diff(dom.sigma).*lhs(:,1);
betS = dom.beta (1) + diff(dom.beta ).*lhs(:,2);
rhoS = dom.rho  (1) + diff(dom.rho  ).*lhs(:,3);
pars = [sigS,betS,rhoS];
lyapSweep = zeros(nSamples,1,'single');
parfor s = 1:nSamples
    lyapSweep(s) = runTraj(pars(s,:), tspanSweep, dtSkip);
end
figure('Name','Lyap map','Color','w');
scatter3(pars(:,1),pars(:,2),pars(:,3),40,lyapSweep,'filled');
xlabel('\sigma'); ylabel('\beta'); zlabel('\rho'); colorbar;
title('Lyapunov spectrum (LHS)');

%% 5 | RHO BIFURCATION SCAN
rhoVec = linspace(dom.rho(1), dom.rho(2), 200);
maxZ   = zeros(size(rhoVec));
minZ   = zeros(size(rhoVec));
for i = 1:numel(rhoVec)
    z_full = ode45(@(t,x)lorenzRHS(t,x,[16 4/3 rhoVec(i)]), ...
                   tspanSweep, [1 1 1]).y(3,:);
    nPts = numel(z_full);
    if nPts >= 1000
        seg = z_full(end-999:end);
    else
        seg = z_full;
    end
    maxZ(i) = max(seg);
    minZ(i) = min(seg);
end
figure('Name','Bifurcation','Color','w');
plot(rhoVec,maxZ,'k', rhoVec,minZ,'k');
xlabel('\rho'); ylabel('Z extrema');
title('Bifurcation diagram (\rho sweep)'); grid on;

%% 6 | UNCERTAINTY QUANTIFICATION
pNom   = [sigma0 beta0 rho0];
lyapUQ = zeros(nEnsembleIC,1);
parfor k = 1:nEnsembleIC
    lyapUQ(k) = runTraj(pNom, tspanSweep, dtSkip, [1 1 1]+0.05*randn(1,3));
end
ci95 = prctile(bootstrp(bootstrapRuns,@mean,lyapUQ), [2.5 97.5]);
fprintf('UQ: Lyap mean %.4f, 95%% CI [%.4f %.4f]\n', mean(lyapUQ), ci95);

%% 7 | POD / SVD REDUCTION
snap = ode45(@(t,x)lorenzRHS(t,x,pNom), tspanSweep, [1 1 1]).y(:,501:end);
[U,S,~] = svd(snap,'econ');
e = diag(S).^2 / sum(diag(S).^2);
fprintf('POD: first two modes capture %.2f%%\n', 100*sum(e(1:2)));
solRed = ode45(@(t,a)U(:,1:2)\lorenzRHS(t,U(:,1:2)*a,pNom), ...
               tspanSweep, U(:,1:2)\[1;1;1]);
figure('Name','POD compare','Color','w');
plot(solRed.y(1,:),solRed.y(2,:),'r'), hold on;
pr = U(:,1:2)'*snap; plot(pr(1,:),pr(2,:),'k.');
legend('reduced','fullprojection'), axis equal;
title('POD 2-mode');

%% 8 | FINITE-DIFF SENSITIVITY
h     = 1e-3;
baseL = runTraj(pNom, tspanSweep, dtSkip);
grad  = zeros(1,3);
for j = 1:3
    dp = pNom; dp(j) = dp(j) + h;
    grad(j) = (runTraj(dp, tspanSweep, dtSkip) - baseL)/h;
end
fprintf('Finite-diff dLyap/d[σ β ρ] = [%.4f %.4f %.4f]\n', grad);

%% 9 | SAVE ALL
res = 'Lorenz_Elite_Out';
if ~exist(res,'dir'), mkdir(res); end
figs = findall(0,'Type','figure');
for F = figs'
    if isgraphics(F,'figure')
        nm = regexprep(get(F,'Name'),'[^\w]','_');
        saveas(F, fullfile(res,[nm '.png']));
    end
end
save(fullfile(res,'workspace.mat'));
fprintf('\nDONE — outputs saved to "%s".\n',res);

%%========================================================================
%% LOCAL FUNCTIONS
%%========================================================================
function dx = lorenzRHS(~,x,p)
    dx = [p(1)*(x(2)-x(1));
          x(1)*(p(3)-x(3)) - x(2);
          x(1)*x(2)       - p(2)*x(3)];
end

function L = runTraj(p,tspan,cut,ic)
    if nargin<4, ic = [1 1 1]; end
    sol = ode45(@(t,x)lorenzRHS(t,x,p), tspan, ic);
    y   = sol.y'; 
    tt  = linspace(tspan(1),tspan(2),size(y,1))';
    idx = tt > cut;
    L   = lyapExponent(y(idx,1), tt(idx));
end

function L = lyapExponent(x,t)
    m   = 3; tau = 10;
    Y   = embed(x,m,tau); Ny = size(Y,1);
    nn  = zeros(Ny,1);
    for i = 1:Ny
        d      = sum((Y - Y(i,:)).^2,2);
        d(i)   = Inf;
        [~,nn(i)] = min(d);
    end
    div = [];
    for i = 1:Ny
        j = nn(i);
        if j+Ny-i>Ny, break; end
        d0 = norm(Y(i,:) - Y(j,:));
        if d0>0, div(end+1)=log(d0); end %#ok<AGROW>
    end
    if numel(div)>1
        p = polyfit(t(1:numel(div)), div', 1);
        L = p(1);
    else
        L = NaN;
    end
end

function Y = embed(x,m,tau)
    N = numel(x); K = N - (m-1)*tau;
    Y = zeros(K,m);
    for k = 1:m
        Y(:,k) = x(1+(k-1)*tau : K+(k-1)*tau);
    end
end
