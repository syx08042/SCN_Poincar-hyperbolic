function d = hyp_dist(r1, theta1, r2, theta2)
% d_H
    delta_theta = pi - abs(pi - abs(theta1 - theta2));
    cosh_d = cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(delta_theta);
    d = acosh(cosh_d);
end
