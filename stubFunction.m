function [Impact_matrix,Impact_heuristic_matrix,Impact_noBus_matrix,...
    userServiceMatrix_MILP,userServiceMatrix_heuristic] = stubFunction(n_delta, ...
    n_simulation,h,lambda_values_string,w,L,G,D,nn,lv,lambda_values,Q,D0,Adj, ...
    DL0,dout,npred0,ne,nde,edge_list,save_lp_model)

    Impact_matrix = zeros(n_simulation,1);
    Impact_heuristic_matrix = zeros(n_simulation,1);
    Impact_noBus_matrix = zeros(n_simulation,1);
    userServiceMatrix_MILP = zeros(n_simulation,1,2);
    userServiceMatrix_heuristic = zeros(n_simulation,1,2);
    hh = 1;
        while hh<=n_simulation
            % fprintf(['Lambda = ' lambda_values_string{lv} '*Q/A | Cycle '...
            %     num2str(h*n_simulation+hh) ' /' num2str(n_delta*n_simulation) '\n'])
            %% User distributions
            dta = 600; % discretization in seconds for the arrival time at a node % 600
            dtr = 600; % discretization in seconds for the ride time % 600
            T = 10*lcm(dta,dtr); % maximum arrival time at final stop
            aa = 0:dta:T-dta;
            rr = dtr:dtr:T; 
            La = length(aa);
            Lr = length(rr);
            
            % for time travel modeling see
            % MODELLING TRAVEL TIME DISTRIBUTION AND ITS INFLUENCE
            % OVER STOCHASTIC VEHICLE SCHEDULING
            [stops_xy,radius] = get_coordinates(); % max radius = 28 km
            radius = 10000; % Legnaro = centr
            [vx,vy] = voronoi(stops_xy(:,1),stops_xy(:,2));
            lambda = lambda_values(lv)*Q / (pi*radius^2);
            M = 5;
            [xx,yy] = inhomogeneous_poisson_generation(radius,-0.8771*1e4,-0.1668*1e4,lambda,M); % Legnaro = centr
            ntu = length(xx); % total number of users
            delta_max = T/min(D0)-1;
            delta = h*delta_max/(n_delta-1);
            user_distributions = zeros(ntu,4); 
            % id number | which node | a | r
            already_served = [];
            for u = 1:ntu
                min_dist_to_stop = +Inf;
                designated_stop = -1;
                for i = 1:nn+1
                    dist_to_stop = norm([xx(u) yy(u)]-stops_xy(i,:));
                    if dist_to_stop < min_dist_to_stop
                        min_dist_to_stop = dist_to_stop;
                        designated_stop = i;
                    end
                end
                Di0 = D0(designated_stop);
                user_distributions(u,1) = u;
                user_distributions(u,2) = designated_stop;
                user_distributions(u,3) = unifrnd(max(0,T-Di0*(1+delta)),T-Di0);
                user_distributions(u,4) = T-user_distributions(u,3);
            end
            
            % disregard those users who are already close to the destination
            ntu_original = ntu;
            u = 0;
            while u < length(user_distributions(:,1))
                u = u + 1;
                if user_distributions(u,2) == nn+1
                    already_served = [already_served user_distributions(u,1)];
                    user_distributions(u,:) = [];
                    ntu = ntu - 1;
                    u = u - 1;
                end
            end
            % if loading_data
            %     load(filename);
            % end
            %-------------------------------------------------------------------------------------
            % already_served
            % fprintf('User distribution (k | node | a [seconds] | r [seconds]):\n\n')
            % user_distributions
            %-------------------------------------------------------------------------------------
            
            users_per_node = zeros(nn,1);
            for i = 1:nn
                users_per_node(i) = length(find(user_distributions(:,2) == i));
            end
            tot_users_per_node = sum(users_per_node);
            
            %% Impact weights
            alpha = 0;
            beta = 1;
            qI = Adj;
            l0 = 18;
            liu = 4.5;
            max_impact = 0;
            for i = 1:nn+1
                for j = 1:nn+1
                    if qI(i,j) > 0
                        users_ij = 0;
                        if i <= nn
                            users_ij = users_per_node(i);
                            if j <= nn     
                                users_ij = users_ij + users_per_node(j);
                            end
                        else
                            if j <= nn
                                users_ij = users_per_node(j);
                            end
                        end
                        Qij = users_ij/Q;
                        qI(i,j) = L(i,j);
                        if l0 > L(i,j)
                            qI(i,j) = +Inf;
                        else
                            max_impact = max_impact + qI(i,j);
                        end
                    else
                        qI(i,j) = -1;
                    end
                end
            end
            GammaI = -ones(nn,1);
            max_impact = max_impact + alpha*sum(GammaI);
            mI = cell(nn,1);
            min_impact = 0;
            for i = 1:nn
                mI{i} = zeros(users_per_node(i),1);
                for k = 1:users_per_node(i)
                    mI{i}(k) = DL0(i);
                    min_impact = min_impact - mI{i}(k);
                end
            end
            min_impact_noBus = min_impact;
            min_impact = beta*min_impact;
            % if loading_data
            %     load(filename);
            %     GammaI = gamma;
            % end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Categories
            [category,theta,nc,ntc] = get_categories(nn,aa,rr,user_distributions);
            
            % computing max_nil for the upper bound of zil
            max_nil = 0;
            for i = 1:nn
                for l = 1:nc(i)
                    nil = category{i}(l,3);
                    if nil > max_nil
                        max_nil = nil;
                    end
                end
            end
            
            
            %% initialization of a feasible solution through heuristics
            % cycling over last stop times + greedy approach + bang-for-buck strategy 
            score = 0;              % number of picked-up users 
            nvis = +Inf;            % number of visited nodes
            P0 = [];                % solution
            d = zeros(nn+1,1)+Inf;  % departure times
            category_index = [];    % category of users for each node
            dt = gcd(dta,dtr);      % discretization time interval
            for t = dt:dt:T         % cycling over last stop times 
                [score_,nvis_,S_,d_,category_index_] =...
                    shuttle_orienteering(t,nn,G,w,nc,category,users_per_node,aa,rr,L);
                if score_ > score || score_ == score && d_(end) < d(end) ||...   
                        score_ == score && d_(end) == d(end) && nvis_ < nvis 
                    score = score_;
                    nvis = nvis_;
                    P0 = S_;
                    d = d_;
                    category_index = category_index_;
                end
            end
            fprintf('H ')
            if rand < 0.1
                fprintf('\n')
            end
            
            
            %% recap after using heuristics
            pickedup_users = zeros(ntu,1);
            users_in_bus = 0;
            user_distrib_ordered = sortrows(user_distributions,3);
            passengersPerStop = [P0(:,1), zeros(length(P0(:,1)),1)];
            for k = 1:ntu
                i = user_distrib_ordered(k,2);
                i_ = find(P0(:,1)==i);
                if ~isempty(i_)
                    aik = ceil(user_distrib_ordered(k,3)/dta)*dta/60;
                    rik = user_distrib_ordered(k,4)/60;
                    if P0(i_,4) >= aik && P0(end,4)-aik <= rik && users_in_bus < Q ...
                           && passengersPerStop(i_,2) < P0(i_,2)              
                        pickedup_users(k) = 1;
                        users_in_bus = users_in_bus + 1;
                        passengersPerStop(i_,2) = passengersPerStop(i_,2)+1;
                    end
                end
            end
            served_users = [pickedup_users...
                user_distrib_ordered(:,1:2) user_distrib_ordered(:,3:4)/60];
            
            userServiceMatrix_heuristic(hh,1,1) = sum(pickedup_users);
            userServiceMatrix_heuristic(hh,1,2) = length(pickedup_users);
            
            %-------------------------------------------------------------------------------------
            % fprintf('Heuristic approach\n')
            % fprintf('User features (is_picked? | k | node | a [minutes] | r [minutes]):\n\n')
            % disp(served_users)
            %-------------------------------------------------------------------------------------
            
            %% MILP formulation: loading the model
            params.Q = Q; % maximum shuttle capacity
            params.T = T;
            params.ntc = ntc;
            params.nc = nc;
            params.category = category;
            params.theta = theta;
            params.nn = nn;
            params.dout = dout;
            params.npred0 = npred0;
            params.w = w;
            params.ne = ne;
            params.nde = nde;
            params.edge_list = edge_list;
            params.D = D;
            params.aa = aa;
            params.rr = rr;
            params.alpha = alpha;
            params.beta = beta;
            params.q = qI;
            params.Gamma = GammaI;
            params.m = mI;
            params.ntu = ntu;
            params.max_nil = max_nil;
            params.min_impact = min_impact;
            params.max_impact = max_impact;
            params.nu = users_per_node;
            [f,intcon,A,b,Aeq,beq,lb,ub,nvar,ncns_eq,ncns_in,varnames,vtype] =...
                milp_model(params,1);
            model.A = A;
            model.b = b;
            model.Aeq = Aeq;
            model.beq = beq;
            model.f = f;
            model.nvar = nvar;
            model.ncns_eq = ncns_eq;
            model.ncns_in = ncns_in;
            model.varnames = varnames;
            model.vtype = vtype;
            
            
            %% build MPS
            %[Contain OK]=BuildMPS(A, b, Aeq, beq, -f, lb, ub, 'shuttle');
            if save_lp_model
                mdl.A = [A; eye(nvar); -eye(nvar); Aeq];
                mdl.rhs = full([b; ub; -lb; beq]);
                mdl.obj = f;
                mdl.varnames = varnames;
                mdl.modelsense = 'min';
                str = '';
                for hhh = 1:ncns_in+2*nvar
                    str = [str '<'];
                end
                for hhh = 1:ncns_eq
                    str = [str '='];
                end
                mdl.sense = str;
                mdl.vtype = vtype;
                mdl.name = ['shuttle_model'];
                gurobi_write(mdl,[mdl.name '.lp']);
                error('STOP')
            end
            
            
            %% MILP settings
            options = optimoptions('intlinprog');
            options.ConstraintTolerance = 1e-4;     % default: 1e-4
            options.CutGeneration = 'advanced';     % default: 'basic'
            options.CutMaxIterations = 10;          % default: 10
            % options.Display = 'iter';               % default: iter
            options.Display = 'off';
            options.Heuristics = 'basic';           % default: 'basic'
            options.HeuristicsMaxNodes = 50;        % default: 50
            options.IntegerPreprocess = 'advanced'; % default: 'basic'
            options.IntegerTolerance = 1e-5;        % default: 1e-5
            options.MaxNodes = 1e7;                 % default: 1e7
            options.MaxTime = 10*60;                % default: 7200
            
            
            %% computation of the MILP solution (no initialization)
            % (i)
            p = intlinprog_grb(f,intcon,A,b,Aeq,beq,lb,ub,[],options);
            Impact = zeros(nvar,1);
            k = 0;
            for i = 1:nn+1
                for j = 1:nn+1
                    if qI(i,j) >= 0
                        k = k + 1;
                        Impact(k) = qI(i,j);
                    end
                end
            end
            for i = 1:nn
                Impact(ne+ntc+i) = alpha*GammaI(i);
            end
            k = ne+ntc+2*nn+1;
            for i = 1:nn
                for u = 1:users_per_node(i)
                    k = k + 1;
                    Impact(k) = -beta*mI{i}(u);
                end
            end
            % (ii)
            params.Impact = Impact'*p;
            stage_1_impact = params.Impact - min_impact;
            [f,intcon,A,b,Aeq,beq,lb,ub,nvar] = milp_model(params,2);
            p = intlinprog_grb(f,intcon,A,b,Aeq,beq,lb,ub,[p; params.Impact],options);
            % (iii)
            Users = zeros(nvar,1);
            k = ne+ntc+2*nn+1;
            for i = 1:nn
                for u = 1:users_per_node(i)
                    k = k + 1;
                    Users(k) = 1;
                end
            end
            params.Impact = [Impact; 0]'*p;
            stage_2_impact = params.Impact - min_impact;
            params.Users = Users'*p;
            [f,intcon,A,b,Aeq,beq,lb,ub,nvar] = milp_model(params,3);
            p = intlinprog_grb(f,intcon,A,b,Aeq,beq,lb,ub,[p; params.Users],options);
            stage_3_impact = [Impact; 0; 0]'*p - min_impact;
            
            %-------------------------------------------------------------------------------------
            % % display solutions
            % fprintf('\nSolution by heuristics (node | #served users | #users at that node | departure [minutes]):\n')
            % disp(P0)
            % score = sum(P0(:,2));
            % fprintf(['#picked-up users through heuristics: ' num2str(score) '\n'])
            %-------------------------------------------------------------------------------------
            %%
            try
    
                [P,v] = MILP2human(p,params);
                
                userServiceMatrix_MILP(hh,1,1) = sum(P(:,2));
                userServiceMatrix_MILP(hh,1,2) = userServiceMatrix_heuristic(hh,1,2);
                %-------------------------------------------------------------------------------------
                % fprintf('\nSolution by MILP solver (node | #served users | #users at that node | departure [minutes]):\n')
                % disp(P)
                % %fprintf('\nSolution by MILP solver (#user id | is picked?):\n')
                % %disp([(1:ntu)' v])
                % score = sum(P(:,2));
                % fprintf(['#picked-up users through NON-initialized MILP: ' num2str(score) '\n'])
                %-------------------------------------------------------------------------------------
                
                %% impact comparison
                % fprintf('\nImpact for the heuristic\n')
                Impact_heuristic = -min_impact;
                for ii = 1:length(P0(:,1))-1
                    i = P0(ii,1);
                    j = P0(ii+1,1);
                    Impact_heuristic = Impact_heuristic + qI(i,j);
                    mI_i = 0;
                    if ~isempty(mI{i})
                        mI_i = mI{i}(1);
                    end
                    Impact_heuristic = Impact_heuristic - beta*passengersPerStop(ii,2)*mI_i;
                end
                %Impact_heuristic
                
                % fprintf('Impact for MILP\n')
                %stage_1_impact   

                Impact_matrix(hh,1) = stage_3_impact;
                Impact_heuristic_matrix(hh,1) = Impact_heuristic;
                Impact_noBus_matrix(hh,1) = -min_impact_noBus;
                hh = hh+1;
            catch ME
                disp(ME.message)
                disp("Errore gestito")
            end
        end % finish cycle with hh
end
