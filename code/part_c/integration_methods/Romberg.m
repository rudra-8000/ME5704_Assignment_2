%%\ Given conductivity data
t = [0 1 2 3 4];
k = [0.10 0.25 0.35 0.40 0.42];

%%\ Romberg integration
n = 3;                 % Number of Romberg levels
R = NaN(n,n);          % Romberg table

% First column: composite trapezoidal rule
for i = 1:n

    % Number of subintervals: 1, 2, 4
    N = 2^(i-1);

    % Select corresponding data points
    idx = 1:(4/N):5;

    % Composite trapezoidal approximation
    R(i,1) = trapz(t(idx), k(idx));

end

% Richardson extrapolation
for j = 2:n
    for i = j:n

        R(i,j) = R(i,j-1) + ...
            (R(i,j-1) - R(i-1,j-1)) / ...
            (4^(j-1) - 1);

    end
end

%%\ Display Romberg table
disp('Romberg Integration Table:'); %[output:9bdbaa66]
disp(R); %[output:9c79f305]

%%\ Final integral
I_romberg = R(n,n);

fprintf('Romberg integral = %.6f\n', I_romberg); %[output:533017da]

%%\ Heat transferred
Qk = 482.7;              % From Part (b)
H = Qk * I_romberg;

fprintf('Heat transferred = %.3f J\n', H); %[output:940a7541]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:9bdbaa66]
%   data: {"dataType":"text","outputData":{"text":"Romberg Integration Table:\n","truncated":false}}
%---
%[output:9c79f305]
%   data: {"dataType":"text","outputData":{"text":"    1.0400       NaN       NaN\n    1.2200    1.2800       NaN\n    1.2600    1.2733    1.2729\n\n","truncated":false}}
%---
%[output:533017da]
%   data: {"dataType":"text","outputData":{"text":"Romberg integral = 1.272889\n","truncated":false}}
%---
%[output:940a7541]
%   data: {"dataType":"text","outputData":{"text":"Heat transferred = 614.423 J\n","truncated":false}}
%---
