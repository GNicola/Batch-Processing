function Amat = AMAT(X);
 % This fuction allow us to compute the A(t) matrix 

 % Dynamical state
 x = X(1);
 y = X(2);
 u = X(3);
 v = X(4);

 % Constants
 mu = X(5);

 % Site location in ECEF
 % X(6:9) = X = [l2 l3 l4 l5]
 
 % Other costants
 Re = X(10);
 % theta = X(12);
 om = X(11);

 % Position magnitude
 r = norm(X(1:2));
 
 % Allocate space for matrix parts
 dV_dR = zeros(2);
 dV_dV = eye(2);

 % Calculate dA_dR
 mufac = (-mu) / r^3;
 defac = 3*mu / r^5;

 du_dx = mufac + defac*x*x;
 du_dy = defac*y*x;

 dv_dx = defac*x*y;
 dv_dy = mufac + defac*y*y;

 dA_dR = [du_dx, du_dy; 
          dv_dx, dv_dy ];

 dA_dV = zeros(2);
 
 % Assemble the matrix
 amat = ...
     [ dV_dR, dV_dV;
       dA_dR, dA_dV ];
 
 Amat = amat;

end


