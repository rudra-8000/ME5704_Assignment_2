%% Validation settings
t = 0:1:4;
dt = t(2) - t(1);
tol = 1e-10;

methods = {'Trapezoidal'; 'Simpson 1/3'; ...
           'Hybrid (Last)'; 'Hybrid (First)'; ...
           'Romberg'; 'Boole'};

% Polynomial degrees used for validation
degrees = [1; 3; 1; 1; 5; 5];

I_exact_all = zeros(6,1);
I_num_all = zeros(6,1);

%% Validate each method
for m = 1:6

    p = degrees(m);

    % Test function f(t) = t^p
    f = t.^p;

    % Exact analytical integral
    I_exact = 4^(p+1)/(p+1);

    switch m

        case 1
            % Composite trapezoidal rule
            I_num = dt/2 * (f(1) + ...
                2*sum(f(2:end-1)) + f(end));

        case 2
            % Composite Simpson's 1/3
            I_num = dt/3 * ...
                (f(1) + 4*f(2) + 2*f(3) ...
                + 4*f(4) + f(5));

        case 3
            % Simpson 3/8 on [0,3]
            % Trapezoidal on [3,4]
            I_38 = 3*dt/8 * ...
                (f(1) + 3*f(2) + 3*f(3) + f(4));

            I_trap = dt/2 * (f(4) + f(5));

            I_num = I_38 + I_trap;

        case 4
            % Trapezoidal on [0,1]
            % Simpson 3/8 on [1,4]
            I_trap = dt/2 * (f(1) + f(2));

            I_38 = 3*dt/8 * ...
                (f(2) + 3*f(3) + 3*f(4) + f(5));

            I_num = I_trap + I_38;

        case 5
            % Three-level Romberg integration
            R = NaN(3,3);

            R(1,1) = trapz(t([1 5]), f([1 5]));
            R(2,1) = trapz(t([1 3 5]), f([1 3 5]));
            R(3,1) = trapz(t, f);

            for j = 2:3
                for i = j:3
                    R(i,j) = R(i,j-1) + ...
                        (R(i,j-1)-R(i-1,j-1)) ...
                        /(4^(j-1)-1);
                end
            end

            I_num = R(3,3);

        case 6
            % Boole's rule
            I_num = 2*dt/45 * ...
                (7*f(1) + 32*f(2) + 12*f(3) ...
                + 32*f(4) + 7*f(5));
    end

    I_exact_all(m) = I_exact;
    I_num_all(m) = I_num;
end

%% Calculate absolute errors
AbsError = abs(I_num_all - I_exact_all);

% Use a scaled tolerance
Passed = AbsError < tol .* max(1,abs(I_exact_all));

%% Display validation results
TestFunction = strcat("t^", string(degrees));

Results = table(methods, TestFunction, ...
    I_exact_all, I_num_all, AbsError, Passed, ...
    'VariableNames', {'Method','TestFunction', ...
    'ExactIntegral','NumericalIntegral', ...
    'AbsoluteError','Passed'});

disp(Results);

%% Overall validation
if all(Passed)
    fprintf('\nAll polynomial validation tests PASSED.\n');
else
    fprintf('\nSome validation tests FAILED.\n');
end