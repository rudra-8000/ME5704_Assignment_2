% A4_TRANSIENT  Reach the steady state by marching dT/dt = A*T - b (alpha = 1) in time.
%   Compares explicit Euler, RK4, implicit Euler and Crank-Nicolson on the cubic-ghost (M3)
%   system; shows the explicit stability limit dt <= 2/|lambda_max| and its violation.
[A, b] = buildSystem('ghost3');  n = numel(b);
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
Tss = gaussElim(A, b);  T0 = 200*ones(n,1);          % plate starts at 200 K, source switched on at t = 0

% ---- extreme eigenvalues by power iteration (A) and inverse iteration (A^-1 via Gauss)
lamMax = powerIter(@(x) A*x, n);                     % largest |lambda| (most negative)
lamMin = 1 / powerIter(@(x) gaussElim(A, x), n);     % smallest |lambda| = slowest decaying mode
dtMax  = 2/abs(lamMax);  tau = 1/abs(lamMin);
fprintf('lambda_max = %.3f, lambda_min = %.3f\n', lamMax, lamMin);
fprintf('explicit Euler limit dt <= %.4f; RK4 limit dt <= %.4f; slowest time constant tau = %.3f\n', dtMax, 2.785/abs(lamMax), tau);

% ---- reference transient: RK4 with a very small step
tEnd = 10*tau;  dtRef = 0.05*dtMax;  nRef = round(tEnd/dtRef);
Tref = marchHeat(A, b, T0, dtRef, nRef, 'rk4');  tRef = (0:nRef)*dtRef;
k22 = 5;                                             % row of node (2,2)  [r = (j-1)*3 + i = 5]

% ---- scheme comparison: steps and accuracy
runs = {'euler', 0.5*dtMax; 'rk4', 0.5*dtMax; 'ieuler', 0.5*dtMax; 'cn', 0.5*dtMax; 'ieuler', 10*dtMax; 'cn', 10*dtMax};
fprintf('%-7s %8s %7s %10s %12s\n', 'scheme', 'dt', 'steps', 'T22(end)', 'max err T22');
H = cell(size(runs,1), 1);  tab = zeros(size(runs,1), 3);   % per run: dt/dtmax, steps, max error
for k = 1:size(runs,1)
    dt = runs{k,2};  ns = round(tEnd/dt);
    H{k} = marchHeat(A, b, T0, dt, ns, runs{k,1});
    tk = (0:ns)*dt;  err = abs(H{k}(k22,:) - interpLin(tRef, Tref(k22,:), tk));
    tab(k,:) = [dt/dtMax, ns, max(err)];
    fprintf('%-7s %8.4f %7d %10.4f %12.4f\n', runs{k,1}, dt, ns, H{k}(k22,end), max(err));
end

% ---- stability limit: explicit Euler just below / just above dtMax
Tstab = zeros(2, 1);  dts = [0.95 1.05]*dtMax;  E = cell(1,2);  tE = cell(1,2);
for k = 1:2
    ns = round(40/dts(k));  Th = marchHeat(A, b, T0, dts(k), ns, 'euler');
    E{k} = max(abs(Th - Tss), [], 1);  tE{k} = (0:ns)*dts(k);
    fprintf('explicit Euler dt = %.2f*dtMax: final max error %.2e\n', dts(k)/dtMax, E{k}(end));
end
% settling time to 0.1 K from the reference run
iS = find(max(abs(Tref - Tss), [], 1) < 0.1, 1);
fprintf('settling time (|T - Tss| < 0.1 K): t = %.2f = %.1f tau\n', tRef(iS), tRef(iS)/tau);
save(fullfile(resDir, 'a4_transient.mat'), 'lamMax', 'lamMin', 'dtMax', 'tau', 'tRef', 'Tref', 'runs', 'tab');

% ---- Figure
f = figure('Visible','off');
subplot(1,2,1);  plot(tRef, Tref(k22,:), 'k-', 'LineWidth', 1.2);  hold on;
sty = {'o', 's', '^', 'd'};  pick = [1 3 5 6];  lab = {'Euler 0.5 dtmax', 'impl. Euler 0.5 dtmax', 'impl. Euler 10 dtmax', 'CN 10 dtmax'};
for q = 1:4
    k = pick(q);  tk = (0:numel(H{k}(k22,:))-1) * runs{k,2};  step = max(1, round(numel(tk)/12));
    plot(tk(1:step:end), H{k}(k22,1:step:end), sty{q}, 'MarkerSize', 4);
end
xlabel('t');  ylabel('$T_{2,2}$ [K]', 'Interpreter', 'latex');  grid on;  set(gca, 'FontSize', 8);
ylim([200 360]);  legend(['RK4 reference', lab], 'Location', 'southeast', 'FontSize', 6);
subplot(1,2,2);  semilogy(tE{1}, E{1}, 'LineWidth', 1.2);  hold on;  semilogy(tE{2}, E{2}, 'LineWidth', 1.2);
xlabel('t');  ylabel('max error to steady state [K]');  grid on;  set(gca, 'FontSize', 8);
legend('0.95 dtmax', '1.05 dtmax', 'Location', 'northwest', 'FontSize', 7);
saveFig(f, 'a_transient', 15, 5.5);

function lam = powerIter(op, n)
% Power iteration: dominant eigenvalue of the linear map OP (Rayleigh quotient).
x = ones(n,1) + 0.1*(1:n)';  x = x/sqrt(x'*x);
for k = 1:500
    y = op(x);  lam = x'*y;  x = y/sqrt(y'*y);
end
end

function v = interpLin(t, y, tq)
% Piecewise-linear interpolation of the samples (t,y) at the query times tq (t uniform).
dt = t(2) - t(1);  k = min(floor(tq/dt) + 1, numel(t)-1);  w = (tq - t(k))/dt;
v = (1-w).*y(k) + w.*y(k+1);
end
