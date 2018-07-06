function Hout = hessianfcn2(x,lambda,ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mu0, lambda0, h)
% Hessian of objective : stiffness matrix
      
        q = x(ddl);%assignment of each degree of freedom to each element
        % Scomposition of the matrix q: in each row there is the list of
        % displacements for each element. The number of rows gives the
        % number of triangles in the mesh.
        q1 = q(:,1);
        q2 = q(:,2);
        q3 = q(:,3);
        q4 = q(:,4);
        q5 = q(:,5);
        q6 = q(:,6);
        q = q.';
        % J is the determinant of the deformation tensor and it's computed
        % for each element.
        J = (b1.*q1+b2.*q3+b3.*q5).*(a1.*q2+a2.*q4+a3.*q6)/4./(AT.^2)-(a1.*q1+a2.*q3+a3.*q5).*(b1.*q2+b2.*q4+b3.*q6)/4./(AT.^2) + (a1.*q2 + a2.*q4 + a3.*q6 + b1.*q1 + b2.*q3 + b3.*q5)/2./AT + 1;        
        H = zeros(max(max(ddl)),max(max(ddl)));
        % This is the loop for assembly
        for i=1:length(a1)
            I = ddl(i,:); % dofs for each triangle
            Kel = 1/J(i)*(lambda0(i)*log(J(i))-mu0(i))*AA(:,:,i)+mu0(i)*B(:,:,i)+1/J(i)^2*(mu0(i)+lambda0(i)*(1-log(J(i))))*(1/AT(i)/2*[b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)].'+AA(:,:,i)*q(:,i))*([b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)]/2/AT(i)+(q(:,i).')*AA(:,:,i));
            H(I,I)= H(I,I)+Kel*h*AT(i); % assembly
        end
        %H=sparse(H);
% Hessian of nonlinear inequality constraint
 %Hg = sparse(size(H));
 %Hout = H + lambda.eqnonlin*Hg;
Hout = H;
