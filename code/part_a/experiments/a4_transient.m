% A4_TRANSIENT  Reach the steady state by marching dT/ds = A*T - b (alpha = 1) in pseudo-time s.
%   Compares explicit Euler, RK4, implicit Euler and Crank-Nicolson on the cubic-ghost (M3)
%   system against the exact solution Tss + exp(A s)(T0 - Tss); shows the explicit
%   stability limit ds <= 2/|lambda_max| and its violation.
[A, b] = buildSystem('ghost3');  n = numel(b);
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
Tss = gaussElim(A, b);  T0 = 200*ones(n,1);          % plate starts at 200 K, source switched on at s = 0

% ---- extreme eigenvalues by power iteration (A) and inverse iteration (A^-1 via Gauss)
lamMax = powerIter(@(x) A*x, n);                     % largest |lambda| (most negative)
lamMin = 1 / powerIter(@(x) gaussElim(A, x), n);     % smallest |lambda| = slowest decaying mode
dtMax  = 2/abs(lamMax);  tau = 1/abs(lamMin);
fprintf('lambda_max = %.3f, lambda_min = %.3f\n', lamMax, lamMin);
fprintf('explicit Euler limit ds <= %.4f; RK4 limit ds <= %.4f; slowest time constant tau = %.3f\n', dtMax, 2.785/abs(lamMax), tau);

% ---- exact reference on a fine grid: one matrix exponential, then repeated multiplication
tEnd = 10*tau;  dtRef = 0.05*dtMax;  nRef = round(tEnd/dtRef);
E1 = expmSeries(A*dtRef);  Tref = zeros(n, nRef+1);  Tref(:,1) = T0;
for k = 1:nRef, Tref(:,k+1) = Tss + E1*(Tref(:,k) - Tss); end
tRef = (0:nRef)*dtRef;
k22 = 5;                                             % row of node (2,2)  [r = (j-1)*3 + i = 5]

% ---- scheme comparison: steps and accuracy (Euler at dt/4 has the same work as RK4 at dt)
runs = {'euler', 0.5*dtMax; 'rk4', 0.5*dtMax; 'ieuler', 0.5*dtMax; 'cn', 0.5*dtMax; 'ieuler', 10*dtMax; 'cn', 10*dtMax; 'euler', 0.125*dtMax};
fprintf('%-7s %8s %7s %10s %12s\n', 'scheme', 'ds', 'steps', 'T22(end)', 'max err T22');
H = cell(size(runs,1), 1);  tab = zeros(size(runs,1), 3);   % per run: ds/dsmax, steps, max error
for k = 1:size(runs,1)
    dt = runs{k,2};  ns = round(tEnd/dt);
    H{k} = marchHeat(A, b, T0, dt, ns, runs{k,1});
    Ex = expmSeries(A*dt);  Tex = T0;  err = 0;          % exact solution on this scheme's own time grid
    for q = 1:ns
        Tex = Tss + Ex*(Tex - Tss);  err = max(err, abs(H{k}(k22,q+1) - Tex(k22)));
    end
    tab(k,:) = [dt/dtMax, ns, err];
    fprintf('%-7s %8.4f %7d %10.4f %12.4f\n', runs{k,1}, dt, ns, H{k}(k22,end), err);
end

% ---- stability limit: explicit Euler just below / just above dtMax
dts = [0.95 1.05]*dtMax;  E = cell(1,2);  tE = cell(1,2);
for k = 1:2
    ns = round(40/dts(k));  Th = marchHeat(A, b, T0, dts(k), ns, 'euler');
    E{k} = max(abs(Th - Tss), [], 1);  tE{k} = (0:ns)*dts(k);
    fprintf('explicit Euler ds = %.2f*dsMax: final max error %.2e\n', dts(k)/dtMax, E{k}(end));
end
% settling time to 0.1 K from the exact solution
iS = find(max(abs(Tref - Tss), [], 1) < 0.1, 1);
fprintf('settling time (|T - Tss| < 0.1 K): s = %.2f = %.1f tau\n', tRef(iS), tRef(iS)/tau);
save(fullfile(resDir, 'a4_transient.mat'), 'lamMax', 'lamMin', 'dtMax', 'tau', 'tRef', 'Tref', 'runs', 'tab');

% ---- Figure
f = figure('Visible','off');
subplot(1,2,1);  plot(tRef, Tref(k22,:), 'k-', 'LineWidth', 1.2);  hold on;
sty = {'o', 's', '^', 'd'};  pick = [1 3 5 6];
lab = {'explicit Euler, $0.5\,\Delta s_{\max}$', 'implicit Euler, $0.5\,\Delta s_{\max}$', 'implicit Euler, $10\,\Delta s_{\max}$', 'Crank-Nicolson, $10\,\Delta s_{\max}$'};
for q = 1:4
    k = pick(q);  tk = (0:numel(H{k}(k22,:))-1) * runs{k,2};
    plot(tk, H{k}(k22,:), ['-' sty{q}], 'MarkerSize', 3.5, 'LineWidth', 0.8);   % marker at every step
end
xlabel('pseudo-time $s$', 'Interpreter', 'latex');  ylabel('$T_{2,2}$ [K]', 'Interpreter', 'latex');  grid on;  set(gca, 'FontSize', 8);
ylim([200 360]);  legend(['exact', lab], 'Interpreter', 'latex', 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 6);   % below the axes, clear of the curves
subplot(1,2,2);  semilogy(tE{1}, E{1}, 'LineWidth', 1.2);  hold on;  semilogy(tE{2}, E{2}, 'LineWidth', 1.2);
xlabel('pseudo-time $s$', 'Interpreter', 'latex');  ylabel('max error to steady state [K]');  grid on;  set(gca, 'FontSize', 8);
legend('$0.95\,\Delta s_{\max}$', '$1.05\,\Delta s_{\max}$', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 7);
saveFig(f, 'a_transient', 15, 7);

function lam = powerIter(op, n)
% Power iteration: dominant eigenvalue of the linear map OP (Rayleigh quotient).
x = ones(n,1) + 0.1*(1:n)';  x = x/sqrt(x'*x);
for k = 1:500
    y = op(x);  lam = x'*y;  x = y/sqrt(y'*y);
end
end
