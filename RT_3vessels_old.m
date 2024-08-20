%% Applying a fractionated RT schedule
% 1D model with 1 vessel inside the domain, 2 at the boundaries and Neumann
% BCs for c
clear all, clc, close all, format compact

%% PDE model
% Import tumour information
tumours = readtable('parameters_1vessel_flux.csv');
q1 = tumours.q_1;
q3 = tumours.q_3;
L_values = tumours.optL;
V0_values = tumours.V_0;
j = 3;

%% Defining the parameters
k = 1e-2;
c_min = 1e-2;
g = 5;
q_1 = q1(j);
q_3 = q3(j);
d_1 = k*q_3;
q_2 = d_1;
L = L_values(j);
d_T = 1e-7*60/L^2;
d_c = 1e-5*60/L^2;
V_0 = V0_values(j);

%% Simulating tumour growth without treatment until convergence for ICs
% Parameters for numerical scheme
dx = V_0/(V_0*1000); % mesh size
x = -1:dx:1; % mesh grid
J = round(2/dx); % number of mesh points
dt = 5000*dx; % time step size
mu = dt/(dx)^2;

% Defining the initial conditions
C0 = ones(J+1,1); T0 = 0.05/(1-V_0)*ones(J+1,1);

%% Find the index j1 and j2 s.t. x(j1)=-V_0 and x(j2)=V_0
j1 = 0; j2 = 0; tol = 1e-15;
for j = 1:J+1
    if abs(x(j) + V_0) < tol
        j1 = j;

    elseif abs(x(j) - V_0) < tol
        j2 = j;
    
    end
end

if j1 == 0 || j2 == 0
    error('No space point at the boundaries of the blood vessel')
end

%% Defining matrices for linear system
A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
LT(j1,j1-1) = 1/(d_T*mu); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
LT(j2,j2-1) = 0; LT(j2,j2) = 0; 
LT(j2,j2+1) = 1/(d_T*mu);
LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;

LC = A; LC(1,1) = -2*(1+dx); LC(1,2) = 2; LC(J+1,J) = 2; LC(J+1,J+1) = -2*(1+dx);
LC(j1-1,j1-1) = 0; LC(j1-1,j1) = -1-2*dx; 
LC(j2+1,j2) = -1-2*dx; LC(j2+1,j2+1) = 0;
LC(j1,j1-1) = (1)/(d_c*mu); LC(j1,j1) = -dx/(d_c*mu); LC(j1,j1+1) = 0;
LC(j2,j2-1) = 0; LC(j2,j2) = -dx/(d_c*mu); LC(j2,j2+1) = (1)/(d_c*mu);

lc = zeros(J+1,1); lc(j1-1) = 2*d_c*mu*dx; lc(j2+1) = 2*d_c*mu*dx; 
lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu*dx; lc(end) = 2*d_c*mu*dx;

I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature

for j = j1+1:j2-1
    LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
    LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
    I1(j,j) = 0; lc(j) = 1; V(j) = 1;
    T0(j) = 0; C0(j) = 1;
end

I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);

L1 = I - d_T*mu*LT; L2 = I - d_c*mu*LC;

% Defining the functions for kinetics
f_T = @(T,C) kinetics_T(T, C, c_min, k, q_1, q_3, d_1);
f_c = @(T,C) kinetics_c(T, C, c_min, k, q_1, q_3, d_1);

%% Simulate tumour growth until convergence
% Initialise vectors
T_old = T0; C_old = C0;
T = L1\(I1*(T_old + dt*f_T(T_old,C_old))); 
C = L2\(I1*(C_old + dt*f_c(T_old,C_old)) + lc);

% Loop through time
tol = 1e-10; m = 2;
while max(abs(T-T_old)) > tol
    T_old = T; C_old = C;
    T = L1\(I1*(T_old + dt*f_T(T_old,C_old)));
    C = L2\(I1*(C_old + dt*f_c(T_old,C_old)) + lc);
    display(m)
    m = m+1;
end

T_eq = T; C_eq = C; % initial conditions for treatment
M = m-1; % number of time steps
t_final_no_treat = M*dt; % final time
time_span_no_treat = 0:dt:t_final_no_treat; % time steps

%% Check convergence
if max(abs(C - C_old)) > tol
    error('The spatial distribution of oxygen concentration has not converged')
end

if max(abs(T - T_old)) > tol
    error('The spatial distribution of tumour density has not converged')
end

%% Simulating a conventional RT schedule

% RT parameters
    % diffusion coefficients
        d_TS = d_T; d_TR = 0; 
    % consumption and proliferation rates 
        theta_1 = 10; theta_2 = 0.1; q_1s = theta_1*q_1; q_3s = theta_2*q_3; 
    % sub-lethally dagamed cells death rate
        d_1s = k*q_3s; q_2s = d_1s;
    % rates of sublethal (nu) and lethal (lambda, lambda_s) RT damage
        nu = 10; lambda = 1; lambda_s = lambda; 
    % RT damage repair rate
        zeta = 5e-3;
    % rate of mitotic catastrophe
        xi = 5e-4;
    % rate of dead cell clearance from tumour
        eta = 5e-5;
    % Number of RT fractions that correspond to the RT doses 
    % 0, 1, 2, 3, 4, 5 Gy, for a total dose of 80 Gy
        num_frac_RT = 20; % [0,80,40,26,20,16]
    % Time between fractions for 5x, 3x, 1x RT per week, respectively
    % 1 day = 1440 min, 2 days = 2880 min, 5 days = 7200 min
        time_btw_frac_RT = 1.44e3; % [1.44e3,2.88e3,7.2e3];
    % Times per week
        times_x_week = 5;
    % Duration of irradiation 10 min
        time_treat = 10;
    % dimensionless RT dose rates that correspond to the RT doses 
    % 0, 1, 2, 3, 4, 5 Gy with R_max = 1
        R = 4/(time_treat); % [0,1,2,3,4,5]/time_treat

% Defining the initial conditions and vasculature
C0 = C_eq; T0 = T_eq; TS0 = zeros(J+1,1); TR0 = zeros(J+1,1);

% Defining the functions for kinetics
h_T = @(T,TS,TR,C,r) RT_kinetics_T(T, TS, TR, C, c_min, V, k, q_1, q_3, g, d_1, ...
                    lambda, nu, zeta, d_1s, lambda_s, xi, eta, q_1s, q_3s, r);
h_TS = @(T,TS,TR,C,r) RT_kinetics_TS(T, TS, TR, C, c_min, V, k, q_1, q_3, g, d_1, ...
                    lambda, nu, zeta, d_1s, lambda_s, xi, eta, q_1s, q_3s, r);
h_TR = @(T,TS,TR,C,r) RT_kinetics_TR(T, TS, TR, C, c_min, V, k, q_1, q_3, g, d_1, ...
                    lambda, nu, zeta, d_1s, lambda_s, xi, eta, q_1s, q_3s, r);
h_c = @(T,TS,TR,C,r) RT_kinetics_c(T, TS, TR, C, c_min, V, k, q_1, q_3, g, d_1, ...
                    lambda, nu, zeta, d_1s, lambda_s, xi, eta, q_1s, q_3s, r);

% Initialise vectors
T = T0; TS = TS0; TR = TR0; C = C0; time_span = 0;

%% Count the RT fraction number
i = 1;
while i < num_frac_RT
    % Irradiation
    % Update time span and paramters for numerical scheme in time
    t_start_irr = time_span(end);
    t_final_irr = t_start_irr + time_treat;
    M_irr = 1000;
    time_span_irr = linspace(t_start_irr,t_final_irr,M_irr);
    dt_irr = time_span_irr(2) - time_span_irr(1);
    mu_irr = dt_irr/(dx^2);
    
    % Redefine matrices for linear system
    A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
    LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
    LT(j1,j1-1) = 1/(d_T*mu_irr); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
    LT(j2,j2-1) = 0; LT(j2,j2) = 0;
    LT(j2,j2+1) = 1/(d_T*mu_irr);
    LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;
    
    LC = A; LC(1,1) = -2*(1+dx); LC(1,2) = 2; LC(J+1,J) = 2; LC(J+1,J+1) = -2*(1+dx); 
    LC(j1-1,j1-1) = 0; LC(j1-1,j1) = -1-2*dx;
    LC(j2+1,j2) = -1-2*dx; LC(j2+1,j2+1) =0;
    LC(j1,j1-1) = (1)/(d_c*mu_irr); LC(j1,j1) = -dx/(d_c*mu_irr); LC(j1,j1+1) = 0;
    LC(j2,j2-1) = 0; LC(j2,j2) = -dx/(d_c*mu_irr); LC(j2,j2+1) = (1)/(d_c*mu_irr);
    
    lc = zeros(J+1,1); lc(j1-1) = 2*d_c*mu_irr*dx; lc(j2+1) = 2*d_c*mu_irr*dx; 
    lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu_irr*dx; lc(end) = 2*d_c*mu_irr*dx;
    
    I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
    V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature
    
    for j = j1+1:j2-1
        LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
        LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
        I1(j,j) = 0; lc(j) = 1; V(j) = 1;
        T0(j) = 0; C0(j) = 1;
    end
    
    I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);
    
    L1_irr = I - d_T*mu_irr*LT; L2_irr = I - d_c*mu_irr*LC;
    L3_irr = I - d_TS*mu_irr*LT; L4_irr = I - d_TR*mu_irr*LT;
    
    % Initial conditions
    T_irr = T(:,end); TS_irr = TS(:,end); TR_irr = TR(:,end); C_irr = C(:,end);
    
    % Solve the PDE model
    for m = 1:M_irr-1
        T_irr(:,m+1) = L1_irr\(I1*(T_irr(:,m) + dt_irr*h_T(T_irr(:,m),TS_irr(:,m), ...
            TR_irr(:,m),C_irr(:,m),R)));
        TS_irr(:,m+1) = L3_irr\(I1*(TS_irr(:,m) + dt_irr*h_TS(T_irr(:,m),TS_irr(:,m), ...
            TR_irr(:,m),C_irr(:,m),R)));
        TR_irr(:,m+1) = L4_irr\(I1*(TR_irr(:,m) + dt_irr*h_TR(T_irr(:,m),TS_irr(:,m), ...
            TR_irr(:,m),C_irr(:,m),R)));
        C_irr(:,m+1) = L2_irr\(I1*(C_irr(:,m) + dt_irr*h_c(T_irr(:,m),TS_irr(:,m), ...
            TR_irr(:,m),C_irr(:,m),R)) + lc);
    end   

    % Determing the break between the this fraction and the next
    if (mod(i,5) == 0)
        time_break = 2.88e3;
        M_btw = 10000;
    else
        time_break = 0;
        M_btw = 5000;
    end
    
    % Tumour growth between fractions
    % Update time span
    t_final_btw = t_final_irr + time_btw_frac_RT - time_treat + time_break;
    time_span_btw = linspace(t_final_irr,t_final_btw,M_btw); 
    dt_btw = time_span_btw(2) - time_span_btw(1);
    mu_btw = dt_btw/(dx^2);
    
    % Redefine matrices for linear system
    A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
    LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
    LT(j1,j1-1) = 1/(d_T*mu_btw); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
    LT(j2,j2-1) = 0; LT(j2,j2) = 0;
    LT(j2,j2+1) = 1/(d_T*mu_btw); 
    LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;
    
    LC = A; LC(1,1) = -2*(1+dx); LC(1,2) = 2; LC(J+1,J) = 2; LC(J+1,J+1) = -2*(1+dx); 
    LC(j1-1,j1) = -1-2*dx; LC(j1-1,j1-1) = 0;
    LC(j2+1,j2) = -1-2*dx; LC(j2+1,j2+1) = 0;
    LC(j1,j1-1) = (1)/(d_c*mu_btw); LC(j1,j1) = -dx/(d_c*mu_btw); LC(j1,j1+1) = 0;
    LC(j2,j2-1) = 0; LC(j2,j2) = -dx/(d_c*mu_btw); LC(j2,j2+1) = (1)/(d_c*mu_btw);
    
    lc = zeros(J+1,1); lc(j1-1) = 2*d_c*mu_btw*dx; lc(j2+1) = 2*d_c*mu_btw*dx; 
    lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu_btw*dx; lc(end) = 2*d_c*mu_btw*dx;
 
    I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
    V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature
    
    for j = j1+1:j2-1
        LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
        LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
        I1(j,j) = 0; lc(j) = 1; V(j) = 1;
        T0(j) = 0; C0(j) = 1;
    end
    
    I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);
    
    L1_btw = I - d_T*mu_btw*LT; L2_btw = I - d_c*mu_btw*LC;
    L3_btw = I - d_TS*mu_btw*LT; L4_btw = I - d_TR*mu_btw*LT;
    
    % Update initial condition
    T_btw = T_irr(:,end); TS_btw = TS_irr(:,end); TR_btw = TR_irr(:,end); 
    C_btw = C_irr(:,end); 
    
    % Solve the PDE model
    for m = 1:M_btw-1
        T_btw(:,m+1) = L1_btw\(I1*(T_btw(:,m) + dt_btw*h_T(T_btw(:,m),TS_btw(:,m), ...
            TR_btw(:,m),C_btw(:,m),0)));
        TS_btw(:,m+1) = L3_btw\(I1*(TS_btw(:,m) + dt_btw*h_TS(T_btw(:,m),TS_btw(:,m), ...
            TR_btw(:,m),C_btw(:,m),0)));
        TR_btw(:,m+1) = L4_btw\(I1*(TR_btw(:,m) + dt_btw*h_TR(T_btw(:,m),TS_btw(:,m), ...
            TR_btw(:,m),C_btw(:,m),0)));
        C_btw(:,m+1) = L2_btw\(I1*(C_btw(:,m) + dt_btw*h_c(T_btw(:,m),TS_btw(:,m), ...
            TR_btw(:,m),C_btw(:,m),0)) + lc);
    end   
    
    % Appending the solutions for the simulated time periods
    time_span = [time_span,time_span_irr, time_span_btw];
    T = [T, T_irr, T_btw];
    TS = [TS, TS_irr, TS_btw];
    TR = [TR, TR_irr, TR_btw];
    C = [C, C_irr, C_btw];

    % Moving to the next fraction
    display(i)
    i = i+1;

end

%% Final fraction
% Irradiation

% Update time span and paramters for numerical scheme in time
t_start_irr = time_span(end);
t_final_irr = t_start_irr + time_treat;
M_irr = 1000;
time_span_irr = linspace(t_start_irr,t_final_irr,M_irr);
dt_irr = time_span_irr(2) - time_span_irr(1);
mu_irr = dt_irr/(dx^2);

% Redefine matrices for linear system
A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
LT(j1,j1-1) = 1/(d_T*mu_irr); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
LT(j2,j2-1) = 0; LT(j2,j2) = 0;
LT(j2,j2+1) = 1/(d_T*mu_irr); 
LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;

LC = A; LC(1,1) = -2*(1+dx); LC(1,2) = 2; LC(J+1,J) = 2; LC(J+1,J+1) = -2*(1+dx);
LC(j1-1,j1) = -1-2*dx; LC(j1-1,j1-1) = 0;
LC(j2+1,j2) = -1-2*dx; LC(j2+1,j2+1) = 0;
LC(j1,j1-1) = (1)/(d_c*mu_irr); LC(j1,j1) = -dx/(d_c*mu_irr); LC(j1,j1+1) = 0;
LC(j2,j2-1) = 0; LC(j2,j2) = -dx/(d_c*mu_irr); LC(j2,j2+1) = (1)/(d_c*mu_irr);

lc = zeros(J+1,1); lc(j1-1) = 2*d_c*mu_irr*dx; lc(j2+1) = 2*d_c*mu_irr*dx; 
lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu_irr*dx; lc(end) = 2*d_c*mu_irr*dx;

I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature

for j = j1+1:j2-1
    LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
    LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
    I1(j,j) = 0; lc(j) = 1; V(j) = 1;
    T0(j) = 0; C0(j) = 1;
end

I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);

L1_irr = I - d_T*mu_irr*LT; L2_irr = I - d_c*mu_irr*LC;
L3_irr = I - d_TS*mu_irr*LT; L4_irr = I - d_TR*mu_irr*LT;

% Initial conditions
T_irr = T(:,end); TS_irr = TS(:,end); TR_irr = TR(:,end); C_irr = C(:,end);

% Solve the PDE model
for m = 1:M_irr-1
    T_irr(:,m+1) = L1_irr\(I1*(T_irr(:,m) + dt_irr*h_T(T_irr(:,m),TS_irr(:,m), ...
        TR_irr(:,m),C_irr(:,m),R)));
    TS_irr(:,m+1) = L3_irr\(I1*(TS_irr(:,m) + dt_irr*h_TS(T_irr(:,m),TS_irr(:,m), ...
        TR_irr(:,m),C_irr(:,m),R)));
    TR_irr(:,m+1) = L4_irr\(I1*(TR_irr(:,m) + dt_irr*h_TR(T_irr(:,m),TS_irr(:,m), ...
        TR_irr(:,m),C_irr(:,m),R)));
    C_irr(:,m+1) = L2_irr\(I1*(C_irr(:,m) + dt_irr*h_c(T_irr(:,m),TS_irr(:,m), ...
        TR_irr(:,m),C_irr(:,m),R)) + lc);
    display(m)
end   

% Post-RT tumour growth
time_break_1 = 2.88e3; 
time_break_2 = 2e5;    

% Update time span
t_final_btw = t_final_irr + time_btw_frac_RT - time_treat + time_break_1;
M_btw = 10000;
time_span_btw = linspace(t_final_irr,t_final_btw,M_btw); 
dt_btw = time_span_btw(2) - time_span_btw(1);
mu_btw = dt_btw/(dx^2);

% Redefine matrices for linear system
A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
LT(j1,j1-1) = 1/(d_T*mu_btw); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
LT(j2,j2-1) = 0; LT(j2,j2) = 0;
LT(j2,j2+1) = 1/(d_T*mu_btw);
LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;

LC = A; LC(1,1) = -2*(1+dx); LC(1,2) = 2; LC(J+1,J) = 2; LC(J+1,J+1) = -2*(1+dx);
LC(j1-1,j1) = -1-2*dx; LC(j1-1,j1-1) = 0;
LC(j2+1,j2) = -1-2*dx; LC(j2+1,j2+1) = 0;
LC(j1,j1-1) = (1)/(d_c*mu_btw); LC(j1,j1) = -dx/(d_c*mu_btw); LC(j1,j1+1) = 0;
LC(j2,j2-1) = 0; LC(j2,j2) = -dx/(d_c*mu_btw); LC(j2,j2+1) = (1)/(d_c*mu_btw);

lc = zeros(J+1,1); lc(j1-1) = 2*d_c*mu_btw*dx; lc(j2+1) = 2*d_c*mu_btw*dx; 
lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu_btw*dx; lc(end) = 2*d_c*mu_btw*dx;

I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature

for j = j1+1:j2-1
    LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
    LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
    I1(j,j) = 0; lc(j) = 1; V(j) = 1;
    T0(j) = 0; C0(j) = 1;
end

I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);

L1_btw = I - d_T*mu_btw*LT; L2_btw = I - d_c*mu_btw*LC;
L3_btw = I - d_TS*mu_btw*LT; L4_btw = I - d_TR*mu_btw*LT;

% Update initial condition
T_btw = T_irr(:,end); TS_btw = TS_irr(:,end); TR_btw = TR_irr(:,end); 
C_btw = C_irr(:,end); 

% Solve the PDE model
for m = 1:M_btw-1
    T_btw(:,m+1) = L1_btw\(I1*(T_btw(:,m) + dt_btw*h_T(T_btw(:,m),TS_btw(:,m), ...
        TR_btw(:,m),C_btw(:,m),0)));
    TS_btw(:,m+1) = L3_btw\(I1*(TS_btw(:,m) + dt_btw*h_TS(T_btw(:,m),TS_btw(:,m), ...
        TR_btw(:,m),C_btw(:,m),0)));
    TR_btw(:,m+1) = L4_btw\(I1*(TR_btw(:,m) + dt_btw*h_TR(T_btw(:,m),TS_btw(:,m), ...
        TR_btw(:,m),C_btw(:,m),0)));
    C_btw(:,m+1) = L2_btw\(I1*(C_btw(:,m) + dt_btw*h_c(T_btw(:,m),TS_btw(:,m), ...
        TR_btw(:,m),C_btw(:,m),0)) + lc);
    display(m)
end   

% Appending the solutions for the simulated time periods
time_span = [time_span,time_span_irr, time_span_btw];
T = [T, T_irr, T_btw];
TS = [TS, TS_irr, TS_btw];
TR = [TR, TR_irr, TR_btw];
C = [C, C_irr, C_btw];

%% Post-treatment tumour growth (comment this out for short-term
%  % response only)
% 
%  % Update time span
%  t_final_growth = time_span(end) + time_break_2;
%  M_growth = 70000;
%  time_span_growth = linspace(time_span(end),t_final_growth,M_growth); 
%  dt_growth = time_span_growth(2) - time_span_growth(2); 
%  mu_growth = dt_growth/(dx^2);
% 
% % Redefine matrices for linear system
% A = -2*eye(J+1)+diag(ones(J,1),-1)+diag(ones(J,1),1);
% LT = A; LT(1,2) = 2; LT(J+1,J) = 2; LT(j1-1,j1) = 0; LT(j1-1,j1-1) = -1;
% LT(j1,j1-1) = 1/(d_T*mu_growth); LT(j1,j1) = 0; LT(j1,j1+1) = 0;
% LT(j2,j2-1) = 0; LT(j2,j2) = 0; % 2/(d_T*mu); 
% LT(j2,j2+1) = 1/(d_T*mu_growth); % -1/(d_T*mu); 
% LT(j2+1,j2) = 0; LT(j2+1,j2+1) = -1;
% 
% LC = A; LC(1,2) = 2*(1-dx); LC(J+1,J) = 2*(1-dx); 
% LC(j1-1,j1) = 0; LC(j1-1,j1-1) = -dx-1;
% LC(j2+1,j2) = 0; LC(j2+1,j2+1) = -dx-1;
% LC(j1,j1-1) = (1-dx)/(d_c*mu_growth); LC(j1,j1) = 0; LC(j1,j1+1) = 0;
% LC(j2,j2-1) = 0; LC(j2,j2) = 0; LC(j2,j2+1) = (1-dx)/(d_c*mu_growth);
% lc = zeros(J+1,1); lc(j1-1) = d_c*mu_growth*dx; lc(j2+1) = d_c*mu_growth*dx; 
% lc(j1) = dx; lc(j2) = dx; lc(1) = 2*d_c*mu_growth*dx; lc(end) = 2*d_c*mu_growth*dx;
% 
% I1 = eye(J+1); I1(j1,j1) = 0; I1(j2,j2) = 0;
% V = zeros(J+1,1); V(j1) = 1; V(j2) = 1; % vasculature
% 
% for j = j1+1:j2-1
%     LC(j,j-1) = 0; LC(j,j) = 0; LC(j,j+1) = 0;
%     LT(j,j-1) = 0; LT(j,j) = 0; LT(j,j+1) = 0;
%     I1(j,j) = 0; lc(j) = 1; V(j) = 1;
%     T0(j) = 0; C0(j) = 1;
% end
% 
% I = sparse(eye(J+1)); I1 = sparse(I1); LC = sparse(LC); LT = sparse(LT);
% 
% L1_growth = I - d_T*mu_growth*LT; L2_growth = I - d_c*mu_growth*LC;
% L3_growth = I - d_TS*mu_growth*LT; L4_growth = I - d_TR*mu_growth*LT;
% 
%  % Update initial conditions
%  T_growth = T(:,end); TS_growth = TS(:,end); TR_growth = TR(:,end);
%  C_growth = C(:,end);
% 
% % Solve the PDE model
%  for m = 1:M_growth
%      T_growth(:,m+1) = L1_growth\(I1*(T_growth(:,m) + dt_growth*h_T(T_growth(:,m), ...
%          TS_growth(:,m),TR_growth(:,m),C_growth(:,m),0)));
%      TS_growth(:,m+1) = L3_growth\(I1*(TS_growth(:,m) + dt_growth*h_TS(T_growth(:,m), ...
%          TS_growth(:,m),TR_growth(:,m),C_growth(:,m),0)));
%      TR_growth(:,m+1) = L4_growth\(I1*(TR_growth(:,m) + dt_growth*h_TR(T_growth(:,m), ...
%          TS_growth(:,m),TR_growth(:,m),C_growth(:,m),0)));
%      C_growth(:,m+1) = L2_growth\(I1*(C_growth(:,m) + dt_growth*h_c(T_growth(:,m), ...
%          TS_growth(:,m),TR_growth(:,m),C_growth(:,m),0)) + lc);
%      display(m)
%  end   
% 
%  % Appending the solutions for the simulated time periods
%  time_span = [time_span,time_span_growth];
%  T = [T, T_growth];
%  TS = [TS, TS_growth];
%  TR = [TR, TR_growth];
%  C = [C, C_growth];

%% Take average over space
% T_total = trapz(x,T,1)/2;
% TS_total = trapz(x,TS,1)/2;
% TR_total = trapz(x,TR,1)/2;
% C_total = trapz(x,C,1)/2;
x_1 = x(1:j1); T_1 = T(1:j1,:);  C_1 = C(1:j1,:); TS_1 = TS(1:j1,:); TR_1 = TR(1:j1,:);
x_2 = x(j1+1:j2-1); T_2 = T(j1+1:j2-1,:); C_2 = C(j1+1:j2-1,:); TS_2 = TS(j1+1:j2-1,:); TR_2 = TR(j1+1:j2-1,:);
x_3 = x(j2:end); T_3 = T(j2:end,:); C_3 = C(j2:end,:); TS_3 = TS(j2:end,:); TR_3 = TR(j2:end,:);
T_total_1 = trapz(x_1,T_1,1); T_total_2 = trapz(x_2,T_2,1); T_total_3 = trapz(x_3,T_3,1);
TS_total_1 = trapz(x_1,TS_1,1); TS_total_2 = trapz(x_2,TS_2,1); TS_total_3 = trapz(x_3,TS_3,1);
TR_total_1 = trapz(x_1,TR_1,1); TR_total_2 = trapz(x_2,TR_2,1); TR_total_3 = trapz(x_3,TR_3,1);
C_total_1 = trapz(x_1,C_1,1); C_total_2 = trapz(x_2,C_2,1); C_total_3 = trapz(x_3,C_3,1);
T_total = (T_total_1 + T_total_2 + T_total_3)/2;
TS_total = (TS_total_1 + TS_total_2 + TS_total_3)/2;
TR_total = (TR_total_1 + TR_total_2 + TR_total_3)/2;
C_total = (C_total_1 + C_total_2 + C_total_3)/2;
            
%% Plot
figure(1)
set(gca,'FontSize',16,'FontName','Helvetica')
hold on
plot(time_span(105001:end), (T_total(105001:end)+TS_total(105001:end))/T_total(1),'color',"#00BFC4", 'LineWidth', 3)
plot(time_span(105001:end), TR_total(105001:end)/T_total(1),'color',"#C77CFF", 'LineWidth', 3)
xlabel('time, t','FontSize',20,'FontName','Helvetica')
title({['q_1=' num2str(q_1) ', q_3=' num2str(q_3) ', V_0=' num2str(V_0) ...
    ', L=' num2str(L)], [num2str(R*time_treat) 'Gy ' num2str(times_x_week) ...
    'x, time treat=' num2str(time_treat), ', # fractions=' num2str(num_frac_RT) ...
    ', tot dose=' num2str(R*time_treat*num_frac_RT)]},'fontsize',14)

figure(2)
set(gca,'FontSize',16,'FontName','Helvetica')
hold on
plot(time_span(105001:end), (C_total(105001:end)),'color',"#F8766D", 'LineWidth', 3)
xlabel('time, t','FontSize',20,'FontName','Helvetica')
ylabel('(c)','FontSize',20,'FontName','Helvetica')
title({['q_1=' num2str(q_1) ', q_3=' num2str(q_3) ', V_0=' num2str(V_0) ...
    ', L=' num2str(L)], [num2str(R*time_treat) 'Gy ' num2str(times_x_week) ...
    'x, time treat=' num2str(time_treat), ', # fractions=' num2str(num_frac_RT) ...
    ', tot dose=' num2str(R*time_treat*num_frac_RT)]},'fontsize',14)

%% ODE Model
% Define the initial conditions
    tt   = 0;
    TT_2 = 0;
    TT_3 = 0;

    % Space limited
    TT_1 = 1-V_0;
    TT_4 = V_0/(V_0 + (q_1/g)*(1-V_0));

    % % Nutrient limited
    % TT_4 = (c_min*(q_1 -3*q_3+(g/c_min+q_3)*V_0 + sqrt((q_1-3*q_3+ ...
    %     (g/c_min+q_3)*V_0)^2 + 4*q_3*(2*(q_1-q_3)+(g-q_1+q_3)*V_0))))...
    %     /(2*(2*(q_1-q_3)+(g-q_1+q_3)*V_0));
    % TT_1 = 2-V_0-c_min/TT_4;

%% Count the RT fraction number
i = 1;
while i < num_frac_RT
    % Irradiation
    % Update time span
    t_1 = linspace(tt(end),tt(end) + time_treat,20); 
    
    % initial conditions
    IC = [TT_1(end), TT_2(end), TT_3(end), TT_4(end)];             
    
    % Solver options
    options = odeset('AbsTol',1e-10, 'RelTol',1e-10);
    
    % Solve the ODE model
    [tt_1,TT] = ode15s(@(t,T) RT(t, T, V_0, q_2, q_1, q_3, g, q_2s, q_1s, q_3s,...
                                 nu, lambda, lambda_s, zeta, eta, xi,...
                                 d_1, d_1s, c_min, R), t_1, IC, options);

    % Determing the break between the this fraction and the next
    if (mod(i,5) == 0)
        time_break = 2.88e3;
    else
        time_break = 0;
    end

    % Tumour growth between fractions
    % Update time span
    t_2 = linspace(tt_1(end), ...
        tt_1(end) + time_btw_frac_RT - time_treat + time_break, 500); 
    
    % Update initial condition
    IC_2 = [TT(end,1), TT(end,2), TT(end,3),TT(end,4)];  
    
    % Solve the model
    [tt_2,TTb] = ode45(@(t,T) RT(t, T, V_0, q_2, q_1, q_3, g, q_2s, q_1s, q_3s,...
                                 nu, lambda, lambda_s, zeta, eta, xi,...
                                 d_1, d_1s, c_min, 0), t_2, IC_2, options); 
    
    % Appending the solutions for the simulated time periods
    tt = vertcat(tt,[tt_1;tt_2]);
    TT_1 = vertcat(TT_1, [TT(:,1);TTb(:,1)]);
    TT_2 = vertcat(TT_2, [TT(:,2);TTb(:,2)]);
    TT_3 = vertcat(TT_3, [TT(:,3);TTb(:,3)]);
    TT_4 = vertcat(TT_4, [TT(:,4);TTb(:,4)]);

    % Moving to the next fraction
    display(i)
    i = i+1;

end

%% Final fraction
% Irradiation

% Update time span
t_1 = linspace(tt(end),tt(end) + time_treat, 100);  

% Update initial conditions
IC = [TT_1(end), TT_2(end), TT_3(end),TT_4(end)];   

% Solve the model
[tt_1,TT] = ode15s(@(t,T) RT(t, T, V_0, q_2, q_1, q_3, g, q_2s, q_1s, q_3s,...
                             nu, lambda, lambda_s, zeta, eta, xi,...
                             d_1, d_1s, c_min, R), t_1, IC, options);

% Post-RT tumour growth
time_break_1 = 2.88e3; 
time_break_2 = 2e5;    

% Update time span
t_2 = linspace(tt_1(end),...
    tt_1(end) + time_btw_frac_RT - time_treat + time_break_1, 5000); 

% Update initial conditions
IC_2 = [TT(end,1), TT(end,2), TT(end,3),TT(end,4)];                                        

% Solve the ODE model
[tt_2,TTb] = ode45(@(t,T) RT(t, T, V_0, q_2, q_1, q_3, g, q_2s, q_1s, q_3s,...
                             nu, lambda, lambda_s, zeta, eta, xi,...
                             d_1, d_1s, c_min, 0), t_2, IC_2, options); 

% Appending the simulated time periods
tt = vertcat(tt,[tt_1;tt_2]);
TT_1 = vertcat(TT_1, [TT(:,1);TTb(:,1)]);
TT_2 = vertcat(TT_2, [TT(:,2);TTb(:,2)]);
TT_3 = vertcat(TT_3, [TT(:,3);TTb(:,3)]);
TT_4 = vertcat(TT_4, [TT(:,4);TTb(:,4)]);

%% Post-treatment tumour growth (comment this out for short-term
% % response only)
% 
% % Update time span
% t_3 = linspace(tt(end),tt(end) + time_break_2, 500); 
% 
% % Update initial conditions
% IC_3 = [TT_1(end), TT_2(end), TT_3(end),TT_4(end)];  
% 
% % Solve the ODE model
% [tt_3,TTc] = ode45(@(t,T) RT(t, T, V_0, q_2, q_1, q_3, g, q_2s, q_1s, q_3s,...
%                              nu, lambda, lambda_s, zeta, eta, xi,...
%                              d_1, d_1s, c_min, 0), t_3, IC_3, options); 
% 
% % Appending the simulated time periods
% tt = vertcat(tt,tt_3);
% TT_1 = vertcat(TT_1, TTc(:,1));
% TT_2 = vertcat(TT_2, TTc(:,2));
% TT_3 = vertcat(TT_3, TTc(:,3));
% TT_4 = vertcat(TT_4, TTc(:,4));
   
%% Plotting the numerical solution
figure(1)
hold on
% Plot the viable tumour cell volume
    plot(tt(7801:end), (TT_1(7801:end)+TT_2(7801:end))/TT_1(1),'--m','LineWidth', 3)   
% Plot the dead cell volume
    plot(tt(7801:end), TT_3(7801:end)/TT_1(1),'--b', 'LineWidth', 3)
legend('(T+T_S)/T_0 PDE','T_R/T_0 PDE', '(T+T_S)/T_0 ODE','T_R/T_0 ODE', ...
    'FontSize',16,'FontName','Helvetica','location','best')

% Plot the oxygen concentration
figure(2)
hold on
plot(tt(7801:end), (TT_4(7801:end)),'--y', 'LineWidth', 3)
yline((c_min),'--k','LineWidth',2)  
legend('(c) PDE','(c) ODE','(c_{min})','FontSize',16,'FontName','Helvetica','location','best')

