function [ sigmaVM, dsigmaVM2_dUf, sigmax, sigmay, sigmaxy, sigmaz, dsigmax_dUf, dsigmay_dUf, dsigmaxy_dUf, dsigmaz_dUf ] = stress_fcn( U, ddl, free_dofs, lambda0, mu0, a1, a2, a3, b1, b2, b3, AT)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
q = U(ddl);
q1 = q(:,1);
q2 = q(:,2);
q3 = q(:,3);
q4 = q(:,4);
q5 = q(:,5);
q6 = q(:,6);
J = (b1.*q1+b2.*q3+b3.*q5).*(a1.*q2+a2.*q4+a3.*q6)/4./(AT.^2)-(a1.*q1+a2.*q3+a3.*q5).*(b1.*q2+b2.*q4+b3.*q6)/4./(AT.^2) + (a1.*q2 + a2.*q4 + a3.*q6 + b1.*q1 + b2.*q3 + b3.*q5)/2./AT + 1;        

%sigmax = (lambda0.*log(J) + mu0.*((Nux.*q + 1).^2 + Nuy.^2.*q.^2 - 1))./J;
sigmax = (lambda0.*log(J) + mu0.*(((b1.*q1+b2.*q3+b3.*q5)./2./AT + 1).^2 + ((a1.*q1+a2.*q3+a3.*q5)./2./AT).^2 - 1))./J;

%sigmay = (lambda0.*log(J) + mu0.*((Nvy.*q + 1).^2 + Nvx.^2.*q.^2 - 1))./J;
sigmay = (lambda0.*log(J) + mu0.*(((a1.*q2+a2.*q4+a3.*q6)./2./AT + 1).^2 + ((b1.*q2+b2.*q4+b3.*q6)./2./AT).^2 - 1))./J;

%sigmaxy = (mu0.*(Nvx.*q.*(Nux.*q + 1) + Nuy.*q.*(Nvy.*q + 1)))./J;
sigmaxy = (mu0.*(((b1.*q2+b2.*q4+b3.*q6)./2./AT).*((b1.*q1+b2.*q3+b3.*q5)./2./AT + 1) + (a1.*q1+a2.*q3+a3.*q5)./2./AT.*((a1.*q2+a2.*q4+a3.*q6)./2./AT+1)))./J;

sigmaz = (lambda0.*log(J))./J;

sigmaVM = sqrt(.5.*((sigmax-sigmay).^2 + (sigmay-sigmaz).^2 + (sigmaz-sigmax).^2 + 6.*sigmaxy.^2));
%sigmaVM = (.5.*((sigmax-sigmay).^2 + (sigmay-sigmaz).^2 + (sigmaz-sigmax).^2 + 6.*sigmaxy.^2));

desigmax_deq = sparse(length(a1),6);
desigmay_deq = sparse(length(a1),6);
dJ_dUf = sparse(length(a1),6);

dJ_dUf(:,1) = b1./(2.*AT) - (q4.*(a1.*b2 - a2.*b1))./(4.*AT.^2) - (q6.*(a1.*b3 - a3.*b1))./(4.*AT.^2);
dJ_dUf(:,2) = (q3.*(a1.*b2 - a2.*b1) + q5.*(a1.*b3 - a3.*b1))./(4.*AT.^2) + a1./(2.*AT);
dJ_dUf(:,3) = b2./(2.*AT) + (q2.*(a1.*b2 - a2.*b1))./(4.*AT.^2) - (q6.*(a2.*b3 - a3.*b2))./(4.*AT.^2);
dJ_dUf(:,4) = a2./(2.*AT) - (q1.*(a1.*b2 - a2.*b1) - q5.*(a2.*b3 - a3.*b2))./(4.*AT.^2);
dJ_dUf(:,5) = b3./(2.*AT) + (q2.*(a1.*b3 - a3.*b1))./(4.*AT.^2) + (q4.*(a2.*b3 - a3.*b2))./(4.*AT.^2);
dJ_dUf(:,6) = a3./(2.*AT) - (q1.*(a1.*b3 - a3.*b1) + q3.*(a2.*b3 - a3.*b2))./(4.*AT.^2);
dsigmax_dJ = (lambda0-sigmax.*J)./J.^2;
dsigmay_dJ = (lambda0-sigmay.*J)./J.^2;
dsigmaz_dJ = lambda0./J.^2.*(1-log(J));
dsigmaxy_dJ = -sigmaxy./J;

dsigmax_dJ = repmat(dsigmax_dJ,1,6);
dsigmay_dJ = repmat(dsigmay_dJ,1,6);
dsigmaxy_dJ = repmat(dsigmaxy_dJ,1,6);
dsigmaz_dJ = repmat(dsigmaz_dJ,1,6);

desigmax_deq(:,1) = (4.*AT.*b1 + 2.*a1.^2.*q1 + 2.*b1.^2.*q1 + 2.*a1.*a2.*q3 + 2.*a1.*a3.*q5 + 2.*b1.*b2.*q3 + 2.*b1.*b3.*q5).*mu0./J./(4.*AT.^2);
desigmax_deq(:,3) = (4.*AT.*b2 + 2.*a2.^2.*q3 + 2.*b2.^2.*q3 + 2.*a1.*a2.*q1 + 2.*a2.*a3.*q5 + 2.*b1.*b2.*q1 + 2.*b2.*b3.*q5).*mu0./J./(4.*AT.^2);
desigmax_deq(:,5) = (4.*AT.*b3 + 2.*a3.^2.*q5 + 2.*b3.^2.*q5 + 2.*a1.*a3.*q1 + 2.*a2.*a3.*q3 + 2.*b1.*b3.*q1 + 2.*b2.*b3.*q3).*mu0./J./(4.*AT.^2);

desigmay_deq(:,2) = (4.*AT.*a1 + 2.*a1.^2.*q2 + 2.*b1.^2.*q2 + 2.*a1.*a2.*q4 + 2.*a1.*a3.*q6 + 2.*b1.*b2.*q4 + 2.*b1.*b3.*q6).*mu0./J./(4.*AT.^2);
desigmay_deq(:,4) = (4.*AT.*a2 + 2.*a2.^2.*q4 + 2.*b2.^2.*q4 + 2.*a1.*a2.*q2 + 2.*a2.*a3.*q6 + 2.*b1.*b2.*q2 + 2.*b2.*b3.*q6).*mu0./J./(4.*AT.^2);
desigmay_deq(:,6) = (4.*AT.*a3 + 2.*a3.^2.*q6 + 2.*b3.^2.*q6 + 2.*a1.*a3.*q2 + 2.*a2.*a3.*q4 + 2.*b1.*b3.*q2 + 2.*b2.*b3.*q4).*mu0./J./(4.*AT.^2);

desigmaxy_deq = [2.*AT.*a1 + a1.^2.*q2 + b1.^2.*q2 + a1.*a2.*q4 + a1.*a3.*q6 + b1.*b2.*q4 + b1.*b3.*q6, 2.*AT.*b1 + a1.^2.*q1 + b1.^2.*q1 + a1.*a2.*q3 + a1.*a3.*q5 + b1.*b2.*q3 + b1.*b3.*q5, 2.*AT.*a2 + a2.^2.*q4 + b2.^2.*q4 + a1.*a2.*q2 + a2.*a3.*q6 + b1.*b2.*q2 + b2.*b3.*q6, 2.*AT.*b2 + a2.^2.*q3 + b2.^2.*q3 + a1.*a2.*q1 + a2.*a3.*q5 + b1.*b2.*q1 + b2.*b3.*q5, 2.*AT.*a3 + a3.^2.*q6 + b3.^2.*q6 + a1.*a3.*q2 + a2.*a3.*q4 + b1.*b3.*q2 + b2.*b3.*q4, 2.*AT.*b3 + a3.^2.*q5 + b3.^2.*q5 + a1.*a3.*q1 + a2.*a3.*q3 + b1.*b3.*q1 + b2.*b3.*q3];
for i=1:6
    desigmaxy_deq(:,i) = desigmaxy_deq(:,i).*mu0./J./(4.*AT.^2);
end

dsigmaVM2_dsigmax = (2*sigmax - sigmay - sigmaz)./sigmaVM/2;
dsigmaVM2_dsigmay = (2*sigmay - sigmax - sigmaz)./sigmaVM/2;
dsigmaVM2_dsigmaz = (2*sigmaz - sigmay - sigmax)./sigmaVM/2;
dsigmaVM2_dsigmaxy = (6*sigmaxy)./sigmaVM/2;

dsigmaVM2_dsigmax = repmat(dsigmaVM2_dsigmax,1,6);
dsigmaVM2_dsigmay = repmat(dsigmaVM2_dsigmay,1,6);
dsigmaVM2_dsigmaz = repmat(dsigmaVM2_dsigmaz,1,6);
dsigmaVM2_dsigmaxy = repmat(dsigmaVM2_dsigmaxy,1,6);

dsigmax_dUf = desigmax_deq + dsigmax_dJ.*dJ_dUf;
dsigmay_dUf = desigmay_deq + dsigmay_dJ.*dJ_dUf;
dsigmaxy_dUf = desigmaxy_deq + dsigmaxy_dJ.*dJ_dUf;
dsigmaz_dUf = dsigmaz_dJ.*dJ_dUf;

dsigmaVM2_dUf = dsigmaVM2_dsigmax.*dsigmax_dUf + dsigmaVM2_dsigmay.*dsigmay_dUf + dsigmaVM2_dsigmaxy.*dsigmaxy_dUf + dsigmaVM2_dsigmaz.*dsigmaz_dUf;

ik = 1:length(a1);
ik = repmat(ik,1,6);

dsigmaVM2_dUf = sparse(ik,ddl,dsigmaVM2_dUf);
dsigmaVM2_dUf = dsigmaVM2_dUf(:,free_dofs);

end