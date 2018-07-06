function [ GKS, dGKS_dUf, dGKS_drho, g_rel ] = GKS_fcn2( U_15, ddl, free_dofs, lambda0, mu0, a1, a2, a3, b1, b2, b3, AT, rho, sigmalim, nt )
%this function computes the stress constraint: we use a unified approach in
% order to approximate the maximum of all the stress constraints through the
% GKS. 
%   We first compute the microscopic stresses and their derivatives, both 
% on free nodal displacements and densities. 
    [ sigmaVM, dsigmaVM_dUf, sigmax, sigmay, sigmaxy, sigmaz, dsigmax_dUf, dsigmay_dUf, dsigmaxy_dUf, dsigmaz_dUf ] = stress_fcn( U_15, ddl, free_dofs, lambda0, mu0, a1, a2, a3, b1, b2, b3, AT);
    P = 4;                          %very important parameter : the bigger is, the more precise will be the GKS (but slower)
    g_rel = (rho.*(sigmaVM/sigmalim-1)).';         % relaxed constraint
    GKS = max(g_rel) + 1/P*log(sum(exp(P*(g_rel-max(g_rel))))) - log(length(g_rel))/P;          % continuous approximation of the maximum
    [I,J] = meshgrid(1:length(free_dofs),1:nt);

    dg_rel_dUf = repmat(rho,1,length(free_dofs)).*(dsigmaVM_dUf/sigmalim);           %%%%% DANGER
    %dGKS_dUf = zeros(length(rho),length(free_dofs));
    dGKS_dUf = exp(P*(g_rel(J)-max(g_rel))).*dg_rel_dUf;                %DANGER
    % If A is a matrix, then sum(A) returns a row vector containing the sum of each column.
    % If A is a vector, then sum(A) returns the sum of the elements.
    dGKS_dUf = sum(dGKS_dUf)/sum(exp(P*(g_rel-max(g_rel))));
    dg_rel_drho = (sigmaVM/sigmalim -1).';
    dGKS_drho = (exp(P*(g_rel-max(g_rel))).*dg_rel_drho/sum(exp(P*(g_rel-max(g_rel))))).';
end