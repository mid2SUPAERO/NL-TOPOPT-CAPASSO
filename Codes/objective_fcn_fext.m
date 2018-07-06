function [ o, do, U_15, dR_drho, K_15 ] = objective_fcn_fext( x, W, penal, ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mu0, lambda0, h, imposed_dofs, ind, free_dofs  )
%The objective function is the compliance due to an external force
% Since the problem is nonlinear, each single part of the calculus of the
% derivatives depends on the nodal displacements. The derivative of the 
% compliance respect to the displacements is the internal force. Then, we 
% use the adjoint sensitivity method adopted in Bhattacharyya's article.
    rho = W*x;              % filtering
    o = sum(rho.*AT)/sum(AT);
    do = AT/sum(AT);   % derivative related to the rho
    mue = mu0 *rho.^penal;           %expression provided by Bhattacharyya
    lambdae = lambda0 *rho.^penal;
    [ U_15, K_15, ~, ~] = nonlinear_tri_fext( ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h, imposed_dofs, ind );
    q = U_15(ddl);  % assignment of each degree of freedom to each element
    % Scomposition of the matrix q: in each row there is the list of
    % displacements for each element. The number of rows gives the
    % number of triangles in the mesh.
    q1 = q(:,1);
    q2 = q(:,2);
    q3 = q(:,3);
    q4 = q(:,4);
    q5 = q(:,5);
    q6 = q(:,6);
    % J is the determinant of the deformation tensor and it's computed
    % for each element.
    J = (b1.*q1+b2.*q3+b3.*q5).*(a1.*q2+a2.*q4+a3.*q6)/4./(AT.^2)-(a1.*q1+a2.*q3+a3.*q5).*(b1.*q2+b2.*q4+b3.*q6)/4./(AT.^2) + (a1.*q2 + a2.*q4 + a3.*q6 + b1.*q1 + b2.*q3 + b3.*q5)/2./AT + 1;
    q=q.';
    m = length(U_15);
    nt = length(a1);         % n is the number of elements (triangles)
    dFint_drho = sparse(m,nt);
    % The force assembly is the same used in gradientfcn2_fext
    for i=1:nt
            I = ddl(i,:);
            Fel=1/J(i)*(lambda0*log(J(i))-mu0)*([b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)]/2/AT(i)+(q(:,i).')*AA(:,:,i))+mu0*([b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)]/2/AT(i)+(q(:,i).')*B(:,:,i));
            dFint_drho(I,i)= (Fel.')*AT(i)*h*rho(i)^(penal-1)*penal;
    end
    % We use the adjoint sensitivity method adopted in Bhattacharyya's article
    dR_drho = dFint_drho;
    dR_drho = dR_drho(free_dofs,:);
    do = W*do;            % density filtering
end


