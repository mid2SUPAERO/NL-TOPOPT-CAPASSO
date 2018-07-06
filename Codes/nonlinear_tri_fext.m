function [ U_15, K_15, Fint_15, phi_15] = nonlinear_tri_fext( ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h, imposed_dofs, ind )
%% Nonlinear system
% We solve the nonlinear force equilibrium system by treating it as a
% minimization problem of the total energy. The gradient indicates the
% forces and the hessian will be the stiffness matrix.
imposed_force_obj = 7*h; %10*h;
step=imposed_force_obj/1;
imposed_force=0;
U = zeros(max(max(ddl)),1);
% The problem will be solved through the application of the fmincon: these
% are the formalization of the constraints in displacements.
Aeq=sparse(max(max(ddl)),1);
beq=Aeq;
Aeq([imposed_dofs])=1;
Aeq=diag(Aeq);
% This is the function handle pointing at the function which builds the
% stiffness matrix: since that fmincon considers the hessian of the
% augmented lagrangian, a lambda is provided as a multiplier for the
% hessian of the nonlinear constraints. But since the constraints are all
% linear, this term disappears.
f_hessian =@(x,lambda) hessianfcn2(x,lambda,ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h);
options = optimoptions('fmincon','MaxIterations',2000,'Display','none',...
    'SpecifyObjectiveGradient',true,'HessianFcn',f_hessian);
% options.StepTolerance=1e-16;
% options.OptimalityTolerance=1e-10;

%% Loop with the force control 
% The solution is founded by adding a small amount of the total force 
% which should be applied
while abs(imposed_force)<abs(imposed_force_obj)
    imposed_force = imposed_force+step;         %we add a small amount of force
    % The expression of the objective function and its gradient change at 
    % each iteration, since they depend on the applied force.
    f_gradient= @(x) gradientfcn2_fext(x,ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h, imposed_force, ind, imposed_dofs);
    U_old=U;
    % The objective function is the total energy Et
    [U, Et] = fmincon(f_gradient,U_old,[],[],Aeq,beq,-20*imposed_force*ones(size(U_old)),20*imposed_force*ones(size(U_old)),[],options);
end
U_15 = U;
[Et, Ftot_15]=f_gradient(U);
phi_15 = Et + imposed_force*U(2*ind-1);
Fint_15 = Ftot_15;
Fint_15(2*ind-1) = Fint_15(2*ind-1) + imposed_force;
K_15 = f_hessian(U,0);
end
