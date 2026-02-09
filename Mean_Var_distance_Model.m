clear;clc;
%% mean var
A = importdata('SCN_2_Aij.mat');
d_E = importdata('SCN_2_d_E_all.mat');
d_H = importdata('SCN_2_d_H_all.mat');
d_E = d_E.*A;d_H = d_H.*A;
d_E_ture = triu(d_E);
nonZeroElements_E = d_E_ture(d_E_ture ~= 0);
d_H_ture = triu(d_H);
nonZeroElements_H = d_H_ture(d_H_ture ~= 0);
nonZeroElements_EE = (nonZeroElements_E - min(nonZeroElements_E)) ./ (max(nonZeroElements_E) - min(nonZeroElements_E));
nonZeroElements_HH = (nonZeroElements_H - min(nonZeroElements_H)) ./ (max(nonZeroElements_H) - min(nonZeroElements_H));
histogram(nonZeroElements_EE);hold on;histogram(nonZeroElements_HH)

mean_E = mean(nonZeroElements_EE);
variance_E = var(nonZeroElements_EE, 1);
mean_H = mean(nonZeroElements_HH);
variance_H = var(nonZeroElements_HH, 1);

%% 
N = 155;   % can be changed
means = 0.1:0.1:0.8;
vars  = 0.01:0.01:0.08;

folder_name = 'E:\data\d_mats_simulateMeanVar_SCN2';
if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

num_pairs = length(means) * length(vars);
pair_list = zeros(num_pairs, 2);
k = 1;
for mi = 1:length(means)
    for vi = 1:length(vars)
        pair_list(k,:) = [means(mi), vars(vi)];
        k = k + 1;
    end
end
L = N*(N-1)/2; 
d_cell = cell(num_pairs, 1);  
info_mean = zeros(num_pairs,1);
info_var  = zeros(num_pairs,1);
info_num  = (1:num_pairs)';
if isempty(gcp('nocreate'))
    parpool;
end
parfor idx = 1:num_pairs
    mu  = pair_list(idx,1);
    var = pair_list(idx,2);
    sigma = sqrt(var);
    d = normrnd(mu, sigma, L, 1);
    d(d < 0) = 0.001;
    d(d > 1) = 1;
    d_cell{idx} = d;
    info_mean(idx) = mu;
    info_var(idx)  = var;
end

for idx = 1:num_pairs
    d = d_cell{idx};  
    D = zeros(N);
    p = 1;
    for i = 1:N-1
        for j = i+1:N
            D(i,j) = d(p);
            D(j,i) = d(p);
            p = p + 1;
        end
    end

    fname = fullfile(folder_name, sprintf('d_%d.mat', idx));
    save(fname, 'D');


    clear D
    d_cell{idx} = [];
end
% save
T = table(info_num, info_mean, info_var, 'VariableNames', {'num','mean','var'});
index_mat_path = fullfile(folder_name, 'index_info.mat');
index_csv_path = fullfile(folder_name, 'index_info.csv');
save(index_mat_path, 'T');
writetable(T, index_csv_path);

%% heatmap matrix
result = importdata('your_result.mat'); 
mean_all = result(:,1);
var_all  = result(:,2);
R_all    = result(:,3);

select_mean = 0.1:0.1:0.8;
select_var  = 0.01:0.01:0.08;

tol = 1e-6;
idx = ismembertol(mean_all, select_mean, tol) & ...
      ismembertol(var_all,  select_var,  tol);

mean_sel = mean_all(idx);
var_sel  = var_all(idx);
R_sel    = R_all(idx);


mean_vals = unique(mean_sel);
var_vals  = unique(var_sel);

M = nan(length(mean_vals), length(var_vals));

for i = 1:length(mean_vals)
    for j = 1:length(var_vals)
        idx2 = (abs(mean_sel - mean_vals(i)) < tol) & ...
               (abs(var_sel  - var_vals(j))  < tol);
        if any(idx2)
            M(i,j) = R_sel(idx2);
        end
    end
end


