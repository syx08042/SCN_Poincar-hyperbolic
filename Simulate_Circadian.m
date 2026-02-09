clear; clc;

ticker = '2'; 
ticker2 = 'H';
filePath1 = fullfile(['SCN_',ticker, '_Aij.mat']);
A = importdata(filePath1);
filePath2 = fullfile(['SCN_',ticker, '_d_', ticker2, '_all.mat']);
dij = importdata(filePath2);
n_cells = length(A);         

alpha = 1;
gamma = 0.5 * ones(1, n_cells);
a     = 1.0 * ones(1, n_cells);
omega = normrnd(2*pi/24, 0.01, [1, n_cells]); 
params = [gamma; a; omega];

x0y0 = rand(1,2*n_cells);
F0 = mean(x0y0(1:n_cells))* ones(1, n_cells);
state0 = [x0y0, F0];                           

tspan = 0:0.01:100000;                    
g = 0.1;                           

[t_out, v_out, R] = Using_Poincare_R_phase(A, g, state0, tspan, params, dij, alpha);

[meanAmp, varAmp] = find_peaks(v_out, t_out, n_cells);




