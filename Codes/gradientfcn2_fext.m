function [f,g] = gradientfcn2_fext(x,ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h, Fext, ind, imposed_dofs)
% Calculate objective f : total energy 

        q = x(ddl); %assignment of each degree of freedom to each element
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
        
        t = 2*ones(size(J)); % tr(C) where C is the right Cauchy tensor
        for i=1:length(a1)
            t(i) = t(i)+ q(i,:)*B(:,:,i)*q(i,:).';
        end
        t = t+2*(b1.*q1 + a1.*q2 + b2.*q3 + a2.*q4 + b3.*q5 + a3.*q6)/2./AT;
        f=h*(AT.')*(.5*lambdae.*(log(J)).^2 -mue.*log(J)+.5*mue.*(t-2));
        f = f - Fext*x(2*ind-1); % external force applied on the horizontal dof of the node number ind
if nargout > 1 % gradient required : forces
    g=zeros(max(max(ddl)),1);
    q=q.';
    for i=1:length(a1)
            I = ddl(i,:);
            Fel=1/J(i)*(lambdae(i)*log(J(i))-mue(i))*([b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)]/2/AT(i)+(q(:,i).')*AA(:,:,i))+mue(i)*([b1(i) a1(i) b2(i) a2(i) b3(i) a3(i)]/2/AT(i)+(q(:,i).')*B(:,:,i));
            g(I)= g(I)+(Fel.')*AT(i)*h;
    end
    g(2*ind-1) = g(2*ind-1) - Fext;
    g(imposed_dofs)= 0;
end
end
