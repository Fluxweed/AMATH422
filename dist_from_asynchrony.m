function [dp] = dist_from_asynchrony(delta_phi)

dp = zeros(1, length(delta_phi));
for i = 1:length(delta_phi)
    p = delta_phi(i);
    if p < pi
        dp(i) = pi - p;
    else
        dp(i) = p - pi;
    end
end
dp = pi - dp;