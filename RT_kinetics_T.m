% Defining the kinetics function for tumour cells

function f = RT_kinetics_T(T, TS, TR, c, c_min, V, k, q_1, q_3, g, d_1, ...
    lambda, nu, mu, d_1s, lambda_s, xi, eta, q_1s, q_3s, R)

f = zeros(length(c),1); 

E = T + TS + TR;

for i = 1:length(c)

    if c(i) >= c_min
        f(i) = k*q_3*c(i)*T(i)*(1-E(i)) - (lambda*c(i)*R+nu*c(i)*R)*T(i) ...
                + mu*TS(i); % tumour cells

    elseif c(i) < c_min
        f(i) = k*q_3*c(i)*T(i)*(1-E(i)) - (d_1*(c_min-c(i)) + ...
                lambda*c(i)*R+nu*c(i)*R)*T(i) + mu*TS(i); % tumour cells
    end

end

end