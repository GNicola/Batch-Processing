function [Ht,Ystar] = HtMAT(X,site);
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
 Ystar = G(X,site);
 rho = Ystar;
 
 % Calculate the derivate to compose Htilde
 drho_dR(1,1) = (x - R*(c*C - S*s)) / rho;
 drho_dR(1,2) = (y - R*(s*C + c*S)) / rho;

 drho_dV = zeros(1,2); 

 % Compose the H tilde matric
 Ht = [drho_dR, drho_dV];

end
