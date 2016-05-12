function [K_Q] = KMAT(X,site,Ystar);
 % Selection the Ground Station 
 if (site == 1)
   lambda = X(13);
 elseif (site == 2)
   lambda = X(6);
 elseif (site == 3)
   lambda = X(7);
 elseif (site == 4)
   lambda = X(8);
 else (site == 5)
   lambda = X(9);
 end

 % Dynamical state
 x = X(1);
 y = X(2);
 u = X(3);
 v = X(4); 

 % Costants
 R = X(10);
 theta = X(12);
 c = cos(theta);
 s = sin(theta);
 C = cos(lambda); 
 S = sin(lambda);

 % Get range
 rho = Ystar;
 drho_dRs = (-R*x*(-c*S - C*s) - R*y*(+C*c - S*s)) / rho;

 % Compose the K matrix
 if (site == 1)
   K_Q = [0, 0, 0, 0];
 elseif (site == 2)
   K_Q = [drho_dRs, 0, 0, 0];
 elseif (site == 3)
   K_Q = [0, drho_dRs, 0, 0];
 elseif (site == 4)
   K_Q = [0, 0, drho_dRs, 0];
 else (site == 5)
   K_Q = [0, 0, 0, drho_dRs];
 end

end

