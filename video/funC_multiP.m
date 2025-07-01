function dx = funC_multiP(t, x, P)
global r K A B

p = P(1);
num_pulses=P(end);
pt0 = P(2:2+num_pulses-1);
pt1 = P(2+num_pulses:end-1);

%%% Perturbation logic
b = B;
for i = 1:length(pt0)
    if t > pt0(i) && t < pt1(i)
        b = p + B;
        break;
    end
end
%%%

dx = (r * x * (1 - x / K) - b * x / (x + A));
end
