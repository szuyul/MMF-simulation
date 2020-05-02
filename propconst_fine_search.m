%propconst_fine_search use the given transcendental equation and intial guess of the root, beta_coarse
% find accurate root, beta_fine
%
% output:
% beta_fine is an 1 by N array of the carefully optimized roots, which are the propagation constants
%
% input:
% transcend_func_abs is the function whose minimum is to be optimized
% beta_coarse is an 1 by N array of the roughly found propagation constants as a stating point
%
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ beta_fine ] = propconst_fine_search( transcend_func_abs, beta_coarse, LB, UB )

num_beta = numel( beta_coarse );
beta_fine = zeros(size(beta_coarse));

options = optimset('MaxFunEvals', 2e6, 'MaxIter', 1e4, 'TolFun', 1e-16);
for ii = 1:num_beta     
    beta_fine(ii) = fminsearchbnd( transcend_func_abs, beta_coarse(ii), LB(ii), UB(ii), options);
end

end

