function  [dX] = EoM(X);
 % inputs:  X = final state
 % outputs: dX= derivative of final state
 
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
 % Re = X(10);
 % theta = X(12);
 % om = X(11);

 % Position magnitude
 r = norm(X(1:2));
 
 % Innitialize the vector
 dX = zeros(4,1);

 % Integrate velocity to find position 
 dX(1:2) = X(3:4);

 % Integrate acceleration to find velocity
 mufac = (-mu) / r^3;
 dX(3) = mufac*x;
 dX(4) = mufac*y;

end
 

