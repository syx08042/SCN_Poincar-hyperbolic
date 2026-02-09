function [t_out, v_out, R] = Using_Poincare_R_phase(A, g, state0, tspan, params, dij, alpha)

    total_len = length(state0);
    n_cells   = total_len / 3;

    max_step = 0.01;
    discard_time = 20000; 

    options = odeset('MaxStep', 0.01, 'RelTol', 1e-5, 'AbsTol', 1e-6);

    tic;
    [t_all, v_all] = ode45(@(t, state) ...
        Poincare_Model(t, state, params, n_cells, g, dij, A, alpha), ...
        tspan, state0, options);
    toc;

    idx_keep = t_all >= discard_time;

    t_out = t_all(idx_keep);
    v_out = v_all(idx_keep, :);

    R_all = zeros(length(t_out), 1);
    for i = 1:length(t_out)
        x = v_out(i, 1:n_cells);
        y = v_out(i, n_cells+1:2*n_cells);
        theta = atan2(y, x);
        R_all(i) = abs(mean(exp(1i * theta)));
    end
    R = mean(R_all);

end
