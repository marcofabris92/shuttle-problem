# shuttle-problem
Code reproducing the numerical simulations of the paper "Optimization Methods for a Shuttle Bus Routing Problem arising in Urban Mobility"

To change scenario, open stubFunction.m and go to line 28. Then pass parameters on your choice to inhomogeneous_poisson_generation(). E.g.:

% Scenario a: 
radius = 10000; % Legnaro = centr
[vx,vy] = voronoi(stops_xy(:,1),stops_xy(:,2));
lambda = lambda_values(lv)\*Q / (pi\*radius^2);
M = 5;
[xx,yy] = inhomogeneous_poisson_generation(radius,-0.8771\*1e4,-0.1668\*1e4,lambda,M); % Legnaro = centr

% Scenario b:
radius = 10000; % Padua station = decentr
[vx,vy] = voronoi(stops_xy(:,1),stops_xy(:,2));
lambda = lambda_values(lv)\*Q / (pi\*radius^2);
M = 3;
[xx,yy] = disomogeneous_poisson_generation(radius,-1.5162\*1e4,0.6144\*1e4,lambda,M); % Padua station = decentr

In addition, parameter beta can be modified at line 85 of stubFunction.m

Finally, in MAIN.m, simulations paramters can be tuned at line 86.
