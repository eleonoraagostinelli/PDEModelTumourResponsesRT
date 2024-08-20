% Defining the kinetics function for oxygen

function f = RT_kinetics_c(T, TS, TR, c, c_min, V, k, q_1, q_3, g, d_1, ...
    lambda, nu, mu, d_1s, lambda_s, xi, eta, q_1s, q_3s, R)

E = T + TS + TR;

f = - q_1*T.*c - q_3*T.*c.*(ones(length(c),1)-E) ...
        - q_1s*TS.*c - q_3s*TS.*c.*(ones(length(c),1)-E); % oxygen 

end