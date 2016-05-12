% Main of the program
% Description: algotithm to work through the batch

 clear all;
 close all;
 clc;
 format long g;

 % Definition
 mu  = 398600.000;			% km^3/s^2
 Re  = 6378.137;			% km
 om  = 7.29115854579367e-5;		% rad/s
 theta = 0.0; 				% rad

 % Load the range (rho)
 data = load('ranges.dat');
 Time = data(:,1);
 Site = data(:,2);
 Observation = data(:,3)./1000;
 l = length(Observation);

 % Initial conditions 
 R0 = [8285.000; 4783.000];
 V0 = [-4.600; 7.900];
 
 % Ground Station position in ECEF
 l1 = rad(8.0);
 l2 = rad(79.8);
 l3 = rad(151.9);
 l4 = rad(224.1);
 l5 = rad(296.2);
 
 lambda = [l1; l2; l3; l4; l5];
 state = [R0; V0; mu; lambda(2:5)];

 % Constant vector
 C = [Re; om; theta; lambda(1)];
 
 % Save the length of the state e constant vector
 n = length(state);
 c = length(C);
 nc= n + c;
 
 % Define STM 
 phi_X = eye(4);
 phi_P = zeros(4,1);
 phibuff = eye(n); 

 % Integration options
 X0 = [state; C; reshape(phi_X,numel(phi_X),1); phi_P];
 t = unique(Time);
 m = length(t);
 [tb,tf] = DIVTIME(t,m);  
 
 % Number of iteration 
 KK = 3; 
 
 %------------------
 % START ITERATION
 %------------------
 for kk = 1 : KK
 
     % Set-up
     L = zeros(n);
     N = zeros(n,1);
     y = zeros(l,1);

     jj = 1;

     %------------------------
     % INTEGRATION OVER TIME
     %------------------------
     options = odeset('RelTol',1e-12,'AbsTol',1e-12);
     [Tb,Xb] = ode45(@CP,tb,X0,options);
     [Tf,Xf] = ode45(@CP,tf,X0,options);
     
     for ii=1 : 145
      contT(ii,1) = Tb(145+1-ii);
      contX(ii,:) = Xb(145+1-ii,:);
     end

     T  = [contT(1:end-1,1); Tf];
     X = [contX(1:end-1,:); Xf];
     X(:,12) = T(:) .* X(:,11); 

     %-----------------------------------
     % CYCLES THROUGHT THE OBSERVATIONS
     %-----------------------------------
     for ii = 1 : l

            % Update the initial condition
            if (Time(ii) == T(jj)) 
               Xi(:,1) = X(jj,:);
               if (jj < length(T))
                  jj = jj + 1;
               end
            end

            % Read observations
            Y = Observation(ii);
            
            % Extract the STM
            phi_X(1:4,1:4) = reshape(Xi(14:29,1),4,4);
            phi_P(1:4,1) = reshape(Xi(30:33,1),4,1);
            
	        % Get the Htilde and K matrixes to compose the H matrix
            [Htilde,Ystar] = HtMAT(Xi(1:nc),Site(ii)); 
            [K_Q] = KMAT(Xi(1:nc),Site(ii),Ystar);
	        H_X = Htilde*phi_X;
	        H_P = Htilde*phi_P;
              
            H = [H_X, H_P, K_Q];
           
     	    % Calculate obseration deviation
            y(ii,1) = Y - Ystar;

	        % Update H
            L = L + (H')*H;
            N = N + H'*y(ii,1);            
                               
     end

     residuals(:,kk) = y;

     % Solve for the final state deviation		
     xhat = L \ N;  % left-division operator which solves: L*xhat = N
     
     % Update the nominal trajectory
     X0(1:n) = X0(1:n) + xhat;
     initial_state_X0(:,kk) = X0; 
     correction_xhat(:,kk) = xhat; 
 end
 pl = [1:1:l];
 save('X0.txt' ,'initial_state_X0', '-ascii', '-double' );
 save('Corretion.txt','correction_xhat', '-ascii', '-double' );
 save('Residuals.txt','residuals', '-ascii', '-double' );

 %---------------
 % PLOT RESIDUALS
 %---------------
 figure;
 hold on;
 grid on;
  plot(pl,residuals(:,1),'k');
  plot(pl,residuals(:,2),'r');
  plot(pl,residuals(:,3),'m');
 hold off;

print -djpg PLOT.jpg

 %---------------
 % SAVE RESULTS 
 %---------------
 

