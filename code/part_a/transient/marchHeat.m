function Th = marchHeat(A, b, T0, dt, nSteps, scheme)
% MARCHHEAT  Pseudo-transient route to the steady state: integrate  dT/dt = A*T - b  (alpha = 1).
%   The steady state of this ODE system is the solution of A*T = b, so marching in time
%   with an ODE/parabolic-PDE scheme gives Laplace's solution as t -> infinity.
%   SCHEME:  'euler'  explicit Euler            T+ = T + dt*f(T)
%            'rk4'    classical Runge-Kutta 4
%            'ieuler' implicit (backward) Euler (I - dt*A) T+ = T - dt*b
%            'cn'     Crank-Nicolson  (I - dt/2*A) T+ = (I + dt/2*A) T - dt*b
%   Returns Th(:,k+1) = temperatures after k steps (Th(:,1) = T0).
n = numel(T0);  I = eye(n);  f = @(T) A*T - b;
Th = zeros(n, nSteps+1);  Th(:,1) = T0;  T = T0;
for k = 1:nSteps
    switch scheme
        case 'euler',  T = T + dt*f(T);
        case 'rk4'
            k1 = f(T);  k2 = f(T + dt/2*k1);  k3 = f(T + dt/2*k2);  k4 = f(T + dt*k3);
            T = T + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        case 'ieuler', T = gaussElim(I - dt*A, T - dt*b);
        case 'cn',     T = gaussElim(I - dt/2*A, (I + dt/2*A)*T - dt*b);
    end
    Th(:,k+1) = T;
end
end
