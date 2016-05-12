function Bmat = BMAT(X);
 % This fuction allow us to compute the B(t) matrix 

 % Dynamical state
 x = X(1);
 y = X(2);
 u = X(3);
 v = X(4);

 % Constants
 % mu = X(5);

 % Site location in ECEF:
 % X(6:9) = XS = [lam2 lam3 lam4 lam5]
 
 % Other costants
 % Re    = X(10);
 % theta = X(12);
 % om    = X(11);

 % Position magnitude
 r = norm(X(1:2));

 % Allocate space for matrix parts
 dV_dmu = zeros(2,1);
 
 % Calculate dA_dR
 fac = (-1) / r^3;
 du_dmu = fac*x;
 dv_dmu = fac*y;

 dA_dmu = [du_dmu; 
           dv_dmu ];


 % Assemble the matrix
 amat = ...
     [ dV_dmu;
       dA_dmu ];
 
 Bmat = amat;

end

