function rho = G(X,site);
 % Dynamical state
 x = X(1);
 y = X(2);
 R = X(10);

 % Select the Ground Station 
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

 % Constants 
 theta = X(12);
 c = cos(theta);
 s = sin(theta);
 C = cos(lambda); 
 S = sin(lambda);

 xs = R*(C*c - S*s);
 ys = R*(S*c + C*s);

 % Compute rho
 rho = sqrt( (x - xs)^2 + (y - ys)^2);
% rho = sqrt( x^2 + y^2 + R^2 - 2*x*R*(C*c + S*s) - 2*y*R*(s*C - c*S) );

end 

