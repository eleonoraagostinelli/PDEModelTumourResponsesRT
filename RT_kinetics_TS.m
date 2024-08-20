% Defining the kinetics function for sub-lethally damaged tumour cells

function f = RT_kinetics_TS(T, TS, TR, c, c_min, V, k, q_1, q_3, g, d_1, ...
    lambda, nu, mu, d_1s, lambda_s, xi, eta, q_1s, q_3s, R)

f = zeros(length(c),1); % [T, TS, TR, c]

E = T + TS + TR;

for i = 1:length(c)

    if c(i) >= c_min
        f(i) = k*q_3s*c(i)*TS(i)*(1-E(i)) - (lambda_s*c(i)*R + mu + xi)*TS(i) ...
                + nu*c(i)*R*T(i); % sub-lethally damaged tumour cells

    elseif c(i) < c_min
        f(i) = k*q_3s*c(i)*TS(i)*(1-E(i)) - (d_1s*(c_min-c(i)) ...
                + lambda_s*c(i)*R + mu + xi)*TS(i) + nu*c(i)*R*T(i); % sub-lethally damaged tumour cells

    end

end


end