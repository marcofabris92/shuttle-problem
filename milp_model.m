function [f,intcon,A,b,Aeq,beq,lb,ub,...
    nvar,ncns_eq,ncns_in,varnames,vtype] =...
    milp_model(p,step)
% optimization variables: 
% [x y s d v I? U?] -> ne ntc nn nn+1 ntu 1? 1?
nvar = p.ne+p.ntc+2*p.nn+p.ntu+1;
varnames = cell(nvar,1);
vars = 0;
for k = 1:p.ne
    vars = vars + 1;
    e = p.edge_list(k,:);
    varnames{vars} = ['x_' num2str(e(1)) '_' num2str(e(2))];
end
for i = 1:p.nn
    for l = 1:p.nc(i)
        vars = vars + 1;
        varnames{vars} = ['y_' num2str(i) '_' num2str(l)];
    end
end
for i = 1:p.nn
    vars = vars + 1;
    varnames{vars} = ['s_' num2str(i)];
end
for i = 1:p.nn+1
    vars = vars + 1;
    varnames{vars} = ['d_' num2str(i)];
end
for i = 1:p.ntu
    vars = vars + 1;
    varnames{vars} = ['v_' num2str(i)];
end
if step >= 2
    nvar = nvar + 1;
    vars = vars + 1;
    varnames{vars} = 'I';
end
if step >= 3
    nvar = nvar + 1;
    vars = vars + 1;
    varnames{vars} = 'U';
end


% objective function selection
switch step
    case 1
        Impact = zeros(nvar,1);
        k = 0;
        for i = 1:p.nn+1
            for j = 1:p.nn+1
                if p.q(i,j) >= 0
                    k = k + 1;
                    Impact(k) = p.q(i,j);
                end
            end
        end
        for i = 1:p.nn
            Impact(p.ne+p.ntc+i) = p.alpha*p.Gamma(i);
        end
        k = p.ne+p.ntc+2*p.nn+1;
        for i = 1:p.nn
            for u = 1:p.nu(i)
                k = k + 1;
                Impact(k) = -p.beta*p.m{i}(u);
            end
        end
        f = Impact;                 % minimizes impact I
    case 2
        Users = zeros(nvar,1);
        k = p.ne+p.ntc+2*p.nn+1;
        for i = 1:p.nn
            for u = 1:p.nu(i)
                k = k + 1;
                Users(k) = 1;
            end
        end
        f = -Users;                 % maximizes the number of users U
    case 3
        f = zeros(nvar,1);
        f(p.ne+p.ntc+2*p.nn+1) = 1; % minimizes d0
end

% [x y s d v I? U?] -> ne ntc nn nn+1 ntu 1? 1?
intcon = [(1:p.ne+p.ntc+p.nn) ...
    (p.ne+p.ntc+2*p.nn+2:p.ne+p.ntc+2*p.nn+p.ntu+1)];
if step >= 3
    intcon = [intcon p.ne+p.ntc+2*p.nn+p.ntu+3];
end
vtype = [repelem('B',p.ne+p.ntc+p.nn,1); repelem('C',p.nn+1,1);
    repelem('B',p.ntu,1)];
if step >= 2
    vtype = [vtype; 'C'];
end
if step >= 3
    vtype = [vtype; 'I'];
end

% [x y s d v I? U?] -> ne ntc nn nn+1 ntu 1? 1?
lb = zeros(nvar,1);
ub = ones(nvar,1);
d_upper = p.ne+p.ntc+p.nn+1:p.ne+p.ntc+2*p.nn+1;
ub(d_upper) = p.T*ub(d_upper);
if step >= 2
    lb(p.ne+p.ntc+2*p.nn+1+p.ntu+1) = p.min_impact;
    ub(p.ne+p.ntc+2*p.nn+1+p.ntu+1) = p.max_impact;
end
if step >= 3
    ub(p.ne+p.ntc+2*p.nn+1+p.ntu+2) = p.ntu;
end

% set of equality constraints
ncns_eq = p.nn+2; % number of equality constraints
Aeq = zeros(ncns_eq,nvar);
beq = zeros(ncns_eq,1);

% set of inequality constraints
ncns_in = 6*p.nn+p.ne+p.nde+p.npred0+1+p.ntu; % number of inequality constraints
if step >= 2
    ncns_in = ncns_in + 1;
end
if step >= 3
    ncns_in = ncns_in + 1;
end
A = zeros(ncns_in,nvar);
b = zeros(ncns_in,1);


% (1): sum(s_i) == 1
Aeq(p.nn+1,p.ne+p.ntc+1:p.ne+p.ntc+p.nn) = ones(1,p.nn);    % sum(s_i)
beq(p.nn+1) = 1;                                            % 1

% (2): sum(x_i0) == 1
for k = 1:p.ne
    if p.edge_list(k,2) == p.nn+1
        Aeq(p.nn+2,k) = 1;                                  % sum(x_i0)
    end
end
beq(p.nn+2) = 1;                                            % 1

% (3): s_i + sum(x_ji) - sum(x_ij) == 0 for each node i
for i = 1:p.nn
    Aeq(i,p.ne+p.ntc+i) = 1;                                % s_i
    for k = 1:p.ne
        if p.edge_list(k,1) == i
            Aeq(i,k) = -1;                                  % - sum(x_ij)
        end
        if p.edge_list(k,2) == i
            Aeq(i,k) = 1;                                   % sum(x_ji)
        end
    end
end

% (4 - valid ineq.): x_ij + x_ji <= 1
h = 0;
for k = 1:p.ne
    i = p.edge_list(k,1);
    j = p.edge_list(k,2);
    if i < j && p.w(i,j) > 0 && p.w(j,i) > 0 && i <= p.nn && j <= p.nn
        h = h + 1;
        kk = 5*p.nn+p.ne+h;
        A(kk,k) = 1;                                        % x_ij
        k_ = 0;
        done = 0;
        while k_ < p.ne && ~done
            k_ = k_ + 1;
            if norm(p.edge_list(k_,:)-[j i]) < eps
                done = 1;
            end
        end
        A(kk,k_) = 1;                                       % x_ji
        b(kk) = 1;                                          % 1
    end
end

% (5): sum(x_ij) <= 1 for each node i
k = 0;
for i = 1:p.nn
    for k_ = 1:p.dout(i)
        k = k + 1;
        A(i,k) = 1;                                         % sum(x_ij)
    end
    b(i) = 1;                                               % 1
end

% (6): (T-D_i0+w_ij)*x_ij + di-dj <= T-D_i0 for each edge (i,j)
for k = 1:p.ne
    kk = p.nn+k;
    i = p.edge_list(k,1);
    j = p.edge_list(k,2);
    A(kk,k) = p.T-p.D(i,p.nn+1)+p.w(i,j);              % (T-D_i0+w_ij)*x_ij
    A(kk,p.ne+p.ntc+p.nn+i) = 1;                       % di
    A(kk,p.ne+p.ntc+p.nn+j) = -1;                      % -dj
    b(kk) = p.T-p.D(i,p.nn+1);                         % T-D_i0
end

% (7): sum((w_ij+D_j0)*x_ij) + di - d0 <= 0 for each node i
for i = 1:p.nn
    ii = p.nn+p.ne+i;
    A(ii,p.ne+p.ntc+p.nn+i) = 1;            % di
    A(ii,p.ne+p.ntc+2*p.nn+1) = -1;         % d0
    for k = 1:p.ne
        if p.edge_list(k,1) == i
            j = p.edge_list(k,2);
            Dj0 = 0;
            if j <= p.nn
                Dj0 = p.D(j,p.nn+1);
            end
            A(ii,k) = p.w(i,j)+Dj0;         % sum((w_ij+D_j0)*x_ij)
        end
    end
end

% (8): -sum(x_ij) + sum(y_il) <= 0 for each node i
k = 0;
h = 0;
for i = 1:p.nn
    ii = 2*p.nn+p.ne+i;
    for k_ = 1:p.dout(i)
        k = k + 1;
        A(ii,k) = -1;           % -sum(x_ij)
    end
    for h_ = 1:p.nc(i)
        h = h + 1;
        A(ii,p.ne+h) = 1;       % sum(y_il)
    end
end

% (9): sum(a_il*y_il) - di <= 0 for each node i
h = 0;
for i = 1:p.nn
    ii = 3*p.nn+p.ne+i;
    for l = 1:p.nc(i)
        h = h + 1;
        ail = p.category{i}(l,1);
        A(ii,p.ne+h) = ail;             % sum(a_il*y_il)
    end
    A(ii,p.ne+p.ntc+p.nn+i) = -1;       % - di
end

% (10): sum((T-a_il-r_il)*y_il) + d0 <= T for each node i
h = 0;
for i = 1:p.nn
    ii = 4*p.nn+p.ne+i;
    for l = 1:p.nc(i)
        h = h + 1;
        ail = p.category{i}(l,1);
        ril = p.category{i}(l,2);
        A(ii,p.ne+h) = p.T-ail-ril;         % sum((T-a_il-r_il)*y_il)
    end
    A(ii,p.ne+p.ntc+2*p.nn+1) = 1;          % d0
    b(ii) = p.T;                            % T
end

% (11 - valid ineq.): sum((D_j0+w_ij-T)*x_ij) + d_i <= 0 for each node i
for i = 1:p.nn
    ii = 5*p.nn+p.ne+p.nde+i;
    A(ii,p.ne+p.ntc+p.nn+i) = 1;            % di
    for k = 1:p.ne
        if p.edge_list(k,1) == i
            j = p.edge_list(k,2);
            Dj0 = 0;
            if j <= p.nn
                Dj0 = p.D(j,p.nn+1);
            end
            A(ii,k) = Dj0+p.w(i,j)-p.T;     % sum((D_j0+w_ij-T)*x_ij)
        end
    end
end

% (12 - valid ineq.): (T-w_i0)*x_i0 - di + d0 <= T 
% for each i such that i is pred(0)
ii = 6*p.nn+p.ne+p.nde;
for k = 1:p.ne
    if p.edge_list(k,2) == p.nn+1
        ii = ii+1;
        i = p.edge_list(k,1);
        A(ii,k) = p.T-p.w(i,p.nn+1);        % (T-w_i0)*x_i0
        A(ii,p.ne+p.ntc+p.nn+i) = -1;       % - di
        A(ii,p.ne+p.ntc+2*p.nn+1) = 1;      % + d0
        b(ii) = p.T;                        % T
    end
end

% (13): sum(sum(v_iu)) <= Q 
ii = 6*p.nn+p.ne+p.nde+p.npred0+1;
k = 0;
for i = 1:p.nn
    for u = 1:p.nu(i)
        k = k + 1;
        A(ii,p.ne+p.ntc+2*p.nn+1+k) = 1;  % sum(sum(v_iu))
    end
end
b(ii) = p.Q;                              % Q

% (14): v_iu <= sum(y_il*theta^i_ul)  for each node, each user
ii = 6*p.nn+p.ne+p.nde+p.npred0+1;
h = 0;
k = 0;
for i = 1:p.nn
    for u = 1:p.ntu
        if ~isempty(p.theta{i,u})
            ii = ii + 1;
            k = k + 1;
            A(ii,p.ne+p.ntc+2*p.nn+1+k) = 1;     % v_iu
            for l = 1:p.nc(i)
                A(ii,p.ne+h+l) = -p.theta{i,u}(l); % -sum(y_il*theta^i_ul)
            end
        end
    end
    h = h + p.nc(i);
end

% (15): I <= I(i) or I <= I(ii)
if step >= 2
    ii = 6*p.nn+p.ne+p.nde+p.npred0+1+p.ntu+1;
    k = 0;
    for i = 1:p.nn+1
        for j = 1:p.nn+1
            if p.q(i,j) >= 0
                k = k + 1;
                A(ii,k) = p.q(i,j);
            end
        end
    end
    for i = 1:p.nn
        A(ii,p.ne+p.ntc+i) = p.alpha*p.Gamma(i);
    end
    k = p.ne+p.ntc+2*p.nn+1;
    for i = 1:p.nn
        for u = 1:p.nu(i)
            k = k + 1;
            A(ii,k) = -p.beta*p.m{i}(u);
        end
    end
    b(ii) = p.Impact;
end

% (16): -U <= -U(ii)
if step >= 3
    ii = 6*p.nn+p.ne+p.nde+p.npred0+1+p.ntu+2;
    k = p.ne+p.ntc+2*p.nn+1;
    for i = 1:p.nn
        for u = 1:p.nu(i)
            k = k + 1;
            A(ii,k) = -1;
        end
    end
    b(ii) = -p.Users;
end


% return constraints
Aeq = sparse(Aeq);
beq = sparse(beq);
A = sparse(A);
b = sparse(b);

end