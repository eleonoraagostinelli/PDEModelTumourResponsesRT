% Defining the kinetics function for lethally damaged tumour cells

function f = RT_kinetics_TR(T, TS, TR, c, c_min, V, k, q_1, q_3, g, d_1, ...
    lambda, nu, mu, d_1s, lambda_s, xi, eta, q_1s, q_3s, R)

% lethally damaged tumour cells
f = lambda.*c.*R.*T + xi*TS + lambda_s.*c.*R.*TS - eta.*TR;

end