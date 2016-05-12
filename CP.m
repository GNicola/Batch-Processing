function [dX] = CP(t,X0); 
 % Cauchy Problem
 X = X0(1:13);
 
 % Extract current STM
 phi_X = reshape(X0(14:29),4,4);
 phi_P = reshape(X0(30:33),4,1);

 % Integrate Acceleration and Velocity
 dX = zeros(length(X0),1);
 dX(1:4) = EoM(X);

 % Integrate the STM
 dphi_X = AMAT(X)*phi_X;
 dphi_P = AMAT(X)*phi_P + BMAT(X);

 dphiX = reshape(dphi_X,numel(dphi_X),1);
 dphiP = reshape(dphi_P,numel(dphi_P),1);
 dphi = [dphiX; dphiP]; 

 dX(14:33) = dphi;
 
% dX(12) = X(11);

end
