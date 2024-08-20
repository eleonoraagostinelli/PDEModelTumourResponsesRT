% Defining the kinetics function for tumour cells
function f = kinetics_T (U, c, c_min, k, q_1, q_3, d_1, V)

f = zeros(length(c),1);

for i = 1:length(c)

    if c(i) >= c_min
    f(i) = k*q_3.*c(i).*U(i).*(1 - U(i) - V(i));

    elseif c(i) < c_min
    f(i) = k*q_3.*c(i).*U(i).*(1 - U(i) - V(i)) - d_1*U(i).*(c_min - c(i));

    end

end

end
