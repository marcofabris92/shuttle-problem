function [score,nvis,S,d,category_index] =...
    shuttle_orienteering(t,nn,G,w,nc,category,nu,aa,rr,L)

% Take the nearest next stop
% initialization of a feasible solution through heuristics greedy approach
S = double.empty(3,0);      % solution: [routh #users departure times]
visited = zeros(nn+1,1);    % boolean flags to mark visited nodes
i = nn+1;                   % current node 
d = zeros(nn+1,1);          % departure times
d(i) = t;
travel_time = 0;            % time interval from node i to final node
arrival_time = t;           % time the shuttle arrives at node i 
ui = 0;                     % number of users collected at node i
new_edge_found = 0;
category_index = zeros(nn,1); 
nu = [nu; 0];
n_iteration = 0;
while ~new_edge_found
    n_iteration = n_iteration+1;
    visited(i) = 1;
    S = [S [i ui nu(i)]'];
    Ni = predecessors(G,i);
    predi = length(Ni);
    b4bi = -ones(predi,1);
    y = zeros(predi,2);
    for j_ = 1:predi
        j = Ni(j_);
        dj = arrival_time-w(j,i);
        if ~visited(j) && dj >= 0 && dj <= d(i)
            n_star = 0;
            for l = 1:nc(j)
                ajl = category{j}(l,1);
                waiting_time = dj-ajl;
                if waiting_time >= 0
                    rjl = category{j}(l,2);
                    njl = category{j}(l,3);           
                    % Remember: if ak <= aa(ja) && 
                    %    ak+rk > aa(ja) && rr(jr) <= ak+rk-aa(ja)
                    if njl > n_star &&...
                            waiting_time+w(j,i)+travel_time <= rjl 
                        n_star = njl;
                        y(j_,1) = l;
                        y(j_,2) = n_star;
                    end
                end
            end
            b4bi(j_) = (y(j_,2)/sum(nu));
        end
    end
    
    % epsilon-greedy approach
    if rand<=0.05 % random approach
        possible_index = [];
        for jj=1:length(b4bi)
            if b4bi(jj) >= 0
                possible_index = [possible_index, jj];
            end
        end
        if isempty(possible_index)
            j_ = 1;
        else
            jj = randi(length(possible_index));
            j_ = possible_index(jj);
        end
        value = b4bi(j_);
    else % bang-for-the-buck approach
        [value,j_] = max(b4bi);
    
        % if two node have equal weight, then take the nearest one
        if sum(value==b4bi)>=2
           node_min_dist = j_;
           for h_=1:length(b4bi)
                if value <= b4bi(h_) && L(h_,i) < L(node_min_dist,i)
                    node_min_dist = h_;
                end
           end
           j_ = node_min_dist;
        end
    end


    j = Ni(j_);
    if value == -1
        new_edge_found = 1;
    else
        d(j) = arrival_time-w(j,i);
        arrival_time = d(j);
        waiting_time = 0;
        if y(j_,1) > 0
            arrival_time = category{j}(y(j_,1),1);
            waiting_time = d(j)-arrival_time;
            category_index(j) = y(j_,1);
        end
        travel_time = travel_time + (w(j,i)+waiting_time);
        i = j;
        ui = y(j_,2);
    end
end
nvis = sum(visited);
S = [fliplr(S)' zeros(nvis,1)];
for i_ = 1:nvis
    i = S(i_,1);
    S(i_,4) = d(i)/60;
end
zc = 0;
while zc < nvis && S(zc+1,2) == 0 && S(zc+1,1) ~= nn+1
    zc = zc + 1;
    d(S(zc,1)) = 0;
end
nvis = nvis - zc;
S = S(zc+1:end,:);
% score = 0;
% for i=1:size(S,1)-1
%     score = score + w(S(i,1), S(i+1,1));
% end
score = sum(S(:,2));
end

