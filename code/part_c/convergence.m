%% Exact integral
I_exact = exp(4) - 1;

% Number of subintervals
N = [4 8 16 32 64];

E_trap = zeros(size(N));
E_simp = zeros(size(N));
E_boole = zeros(size(N));

%% Calculate numerical errors
for j = 1:length(N)

    n = N(j);
    t = linspace(0,4,n+1);
    f = exp(t);
    dt = 4/n;

    % Composite trapezoidal
    I_trap = trapz(t,f);

    % Composite Simpson's 1/3
    I_simp = dt/3 * ...
        (f(1) + f(end) ...
        + 4*sum(f(2:2:end-1)) ...
        + 2*sum(f(3:2:end-2)));

    % Composite Boole's rule
    I_boole = 0;

    for i = 1:4:n
        I_boole = I_boole + 2*dt/45 * ...
            (7*f(i) + 32*f(i+1) + ...
            12*f(i+2) + 32*f(i+3) ...
            + 7*f(i+4));
    end

    % Absolute errors
    E_trap(j) = abs(I_trap-I_exact);
    E_simp(j) = abs(I_simp-I_exact);
    E_boole(j) = abs(I_boole-I_exact);
end

%% Calculate observed convergence orders
p_trap = log2(E_trap(1:end-1)./E_trap(2:end));
p_simp = log2(E_simp(1:end-1)./E_simp(2:end));
p_boole = log2(E_boole(1:end-1)./E_boole(2:end));

%% Plot convergence
figure;
hold on;

h = 4./N;

loglog(h, E_trap, 'bo-', 'LineWidth', 1.5, ...
    'DisplayName', 'Trapezoidal');

loglog(h, E_simp, 'rs-', 'LineWidth', 1.5, ...
    'DisplayName', 'Simpson 1/3');

loglog(h, E_boole, 'kd-', 'LineWidth', 1.5, ...
    'DisplayName', 'Boole');

xlabel('Step size, h');
ylabel('Absolute error');
title('Convergence of Numerical Integration Methods');

legend('Location','northwest');
grid on;
box on;
hold off;

% Display observed convergence orders
fprintf('\nConvergence Test Results:\n');

fprintf(' N      dt       Trap       Simpson      Boole\n');

for i = 2:length(N)
    fprintf('%3d   %7.4f   %8.4f   %8.4f   %8.4f\n', ...
        N(i), 4/N(i), ...
        p_trap(i-1), p_simp(i-1), p_boole(i-1));
end

% Exact integral
I_exact = exp(4) - 1;

% Romberg levels
n = 7;
R = NaN(n,n);

% First column: trapezoidal approximations
for i = 1:n

    N = 2^(i-1);
    t = linspace(0,4,N+1);
    f = exp(t);

    R(i,1) = trapz(t,f);
end

% Richardson extrapolation
for j = 2:n
    for i = j:n
        R(i,j) = R(i,j-1) + ...
            (R(i,j-1)-R(i-1,j-1)) ...
            /(4^(j-1)-1);
    end
end

% Absolute errors
E = abs(R - I_exact);

% Observed convergence orders
p = NaN(n,n);

for j = 1:3
    for i = j+1:n
        p(i,j) = log2(E(i-1,j)/E(i,j));
    end
end

% Display convergence results
N = 2.^(0:n-1)';
dt = 4./N;

Results = table(N, dt, p(:,1), p(:,2), p(:,3), ...
    'VariableNames', ...
    {'N','StepSize','Order2','Order4','Order6'});

disp(Results);

%% Exact integral
I_exact = exp(4) - 1;

% Number of subintervals
% N-1 must be divisible by 3
N = [4 7 10 13 16 19 25 31];

E_first = zeros(size(N));
E_last = zeros(size(N));

%% Calculate numerical errors
for j = 1:length(N)

    n = N(j);
    t = linspace(0,4,n+1);
    f = exp(t);
    dt = 4/n;

    %% A. Trapezoidal rule on last interval
    I_38_last = 0;

    % Composite Simpson 3/8 on first n-1 intervals
    for i = 1:3:n-1
        I_38_last = I_38_last + 3*dt/8 * ...
            (f(i) + 3*f(i+1) + ...
            3*f(i+2) + f(i+3));
    end

    % Trapezoidal rule on final interval
    I_trap_last = dt/2 * (f(n) + f(n+1));

    I_last = I_38_last + I_trap_last;

    %% B. Trapezoidal rule on first interval
    I_trap_first = dt/2 * (f(1) + f(2));

    I_38_first = 0;

    % Composite Simpson 3/8 on remaining intervals
    for i = 2:3:n-1
        I_38_first = I_38_first + 3*dt/8 * ...
            (f(i) + 3*f(i+1) + ...
            3*f(i+2) + f(i+3));
    end

    I_first = I_trap_first + I_38_first;

    %% Absolute errors
    E_last(j) = abs(I_last-I_exact);
    E_first(j) = abs(I_first-I_exact);

end

%% Observed convergence orders
% Since N does not double, use the general formula
p_last = log(E_last(1:end-1)./E_last(2:end)) ...
    ./ log(N(2:end)./N(1:end-1));

p_first = log(E_first(1:end-1)./E_first(2:end)) ...
    ./ log(N(2:end)./N(1:end-1));

%% Display results
Results = table(N(2:end)', (4./N(2:end))', ...
    p_last', p_first', ...
    'VariableNames', ...
    {'N','StepSize','TrapezoidalLast', ...
    'TrapezoidalFirst'});

disp(Results);

%% Plot convergence
figure;
hold on;

h = 4./N;

loglog(h, E_last, 'bo-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Trapezoidal last');

loglog(h, E_first, 'rs-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Trapezoidal first');

xlabel('Step size, h');
ylabel('Absolute error');
title('Convergence of Hybrid Simpson 3/8-Trapezoidal Rule');

legend('Location','northwest');
grid on;
box on;
hold off;
