function h = kinetics_c(U, c, c_min, k, q_1, q_3, d_1, V, g)

% g*(1 - c).*V
h = g*(1 - c).*V - q_1*c.*U - q_3*c.*U.*(ones(length(c),1)-U-V);

end