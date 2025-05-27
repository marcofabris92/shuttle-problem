% SHUTTLE PROBLEM
% M. Fabris, A. Zanella, R. Roberti, R. Carli, L. Boscolo Berto
% May. 17, 2024

clear all
close all
clc

lambda_values = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4];
lambda_values_string = {'025', '050', '075', '100', '125', '150', '175',...
    '200', '250', '300', '350', '400'};

if length(lambda_values) ~= length(lambda_values_string)
    error("Length of the lambda values is different");
end

for lv = 1:length(lambda_values)
    fprintf('Time: %s | Lambda = %s*Q/A \n', datestr(now, 'HH:MM:SS'),lambda_values_string{lv});

    clearvars -except lv lambda_values lambda_values_string

    save_lp_model = 0;
    loading_data = 0;
    addpath '/home/marcof/gurobi1202/linux64/matlab/';
    
    % shuttle capacity
    Q = 100;
    
    
    %% Topology
  
    % % spatial distances 
     %[m]
    L = 1000*table2array(readtable('Adjacency_matrix.csv'));
    
    % % temporal distances: (i,j) means time from node i to node j
     %[s]
    w = L/(40/3.6); % average speed = 40 km/h
     
    G = digraph(w);
    GL = digraph(L);
    Adj = adjacency(G);
    nn = numnodes(G);
    % plot(G)
    D = distances(G);   %[s]
    DL = distances(GL); %[m]
    D0 = D(:,end);      %[s]
    DL0 = DL(:,end);    %[m]
    DL0M = max(DL0);    %[m]
    D = D(1:end-1,:);
    dout = outdegree(G);
    npred0 = length(predecessors(G,nn));
    nsucc0 = length(successors(G,nn));
    nn = nn - 1; % subtract 1 to remove the final node
    ne = numedges(G) - nsucc0; % subtract # of successors of final node
    DLshort = zeros(nn,1)+Inf;
    for i = 1:nn
        spi0 = shortestpath(GL,i,nn+1);
        for k = 1:length(spi0)-1
            Lk = L(spi0(k),spi0(k+1));
            if Lk < DLshort(i)
                DLshort(i) = Lk;
            end
        end
    end
    edge_list = zeros(ne,2);
    k = 0;
    for i = 1:nn
        for j = 1:nn+1
            if w(i,j) > 0
                k = k + 1;
                edge_list(k,:) = [i j];
            end
        end
    end
    nde = 0; % number of well-defined edges
    for i = 2:nn
        for j = 1:i-1
            if w(i,j) > 0 && w(j,i) > 0
                nde = nde + 1;
            end
        end
    end
     
    %% Cycle for different values of lambda
    n_delta = 5; %10;  
    n_simulation = 10; %100;
    Impact_matrix = zeros(n_simulation,n_delta);
    Impact_heuristic_matrix = zeros(n_simulation,n_delta);
    Impact_noBus_matrix = zeros(n_simulation,n_delta);
    userServiceMatrix_MILP = zeros(n_simulation, n_delta,2);
    userServiceMatrix_heuristic = zeros(n_simulation, n_delta,2);

    parfor h = 0:n_delta-1
        [I_m, I_h_m, I_b_m, u_M, u_h] = stubFunction(n_delta, ...
            n_simulation,h,lambda_values_string,w,L,G,D,nn,lv, ...
            lambda_values,Q,D0,Adj,DL0,dout,npred0,ne,nde,edge_list,save_lp_model);
        
        Impact_matrix(:,h+1) = I_m(:);
        Impact_heuristic_matrix(:,h+1) = I_h_m(:);
        Impact_noBus_matrix(:,h+1) = I_b_m(:);
        userServiceMatrix_MILP(:,h+1,:) = u_M(:,1,:);
        userServiceMatrix_heuristic(:,h+1,:) = u_h(:,1,:);
    end 
        
    savepath = sprintf('workspace1c/wrkspc_C_%sQ.mat', lambda_values_string{lv});

    save(savepath)

end % finish cycle with lv

disp("Simulations concluded!!")
