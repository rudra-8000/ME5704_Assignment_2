function [Test, se] = walkGrid(A, b, r0, nWalks)
% WALKGRID  Monte Carlo (random-walk) solution of A*T = b at node R0.
%   Row r of the system reads  T_r = sum_k p_rk T_k + q_r * Tbar_r  with
%       p_rk   = A_rk / (-A_rr)        probability of stepping to neighbour k
%       q_r    = 1 - sum_k p_rk        probability of hitting the edge/line instead
%       Tbar_r = (-b_r/(-A_rr)) / q_r  mean edge/line temperature seen from r
%   so a walker that steps with probabilities p and, when it hits the edge/line, scores
%   Tbar_r has expected score T_r0.  Estimate = mean score over NWALKS walkers.
n = size(A,1);  d = -diag(A);
Pm = A ./ d;  Pm(1:n+1:end) = 0;                 % off-diagonal step probabilities
cum = cumsum(Pm, 2);  q = 1 - cum(:,end);  Tbar = (-b ./ d) ./ q;
s = r0 * ones(nWalks, 1);  score = zeros(nWalks, 1);  alive = true(nWalks, 1);
while any(alive)
    u    = rand(nWalks, 1);
    next = 1 + sum(u > cum(s,:), 2);             % categorical draw; n+1 = hit edge/line
    stop = alive & (next > n);
    score(stop) = Tbar(s(stop));                 % absorbed: score the boundary average
    alive = alive & ~stop;
    s(alive) = next(alive);
end
Test = mean(score);  se = std(score) / sqrt(nWalks);
end
