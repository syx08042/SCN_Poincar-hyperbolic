function rates = Poincare_Model(t, statevars, params, n_cells, g, dij, A, alpha)
    x = statevars(1:n_cells)';
    y = statevars(n_cells+1:2*n_cells)';
    F = statevars(2*n_cells+1:3*n_cells)';

    gamma = params(1,:);  
    a = params(2,:);
    omiga = params(3,:); 
    
    [i_idx, j_idx] = find(A);
    weights = (A(sub2ind(size(A), i_idx, j_idx)) ./ dij(sub2ind(size(dij), i_idx, j_idx))).^alpha;
    weights(isinf(weights) | isnan(weights)) = 0;
    edges = [i_idx, j_idx, weights];
    W = (A ./ dij).^alpha;
    W(isinf(W)) = 0; W(isnan(W)) = 0;
    s = sum(W, 1)';  
    s(s == 0) = 1;   
    rates = zeros(3*n_cells,1);
    r = sqrt(x.^2 + y.^2);
    % dx/dt
    rates(1:n_cells) = gamma .* x .* (a - r) - omiga .* y + g .* F;
    % dy/dt
    rates(n_cells+1:2*n_cells) = gamma .* y .* (a - r) + omiga .* x;
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        w_ij = edges(k, 3);
        rates(2*n_cells + i) = rates(2*n_cells + i) + w_ij * rates(j);
    end
    rates(2*n_cells+1:3*n_cells) = rates(2*n_cells+1:3*n_cells)./ s;
    
end