function [P,v] = MILP2human(pp,p)
% reconstructing the routh information
Lp = int64(1+sum(pp(1:p.ne)));                  % path length
P = zeros(Lp,4);                                % MILP solution
i = find(pp(p.ne+p.ntc+1:p.ne+p.ntc+p.nn));     % first node
i_ = 1;
all_nodes_found = 0;
v = pp(p.ne+p.ntc+2*p.nn+2:p.ne+p.ntc+2*p.nn+1+p.ntu);
while ~all_nodes_found
    k = 0;
    new_edge_found = 0;
    while ~new_edge_found
        k = k + 1;
        if int64(pp(k)) == 1 && p.edge_list(k,1) == i
            new_edge_found = 1;
            P(i_,1) = i;
            j = sum(p.nu(1:i-1));
            P(i_,2) = sum(v(j+1:j+p.nu(i)));
            P(i_,3) = p.nu(i);
            P(i_,4) = pp(p.ne+p.ntc+p.nn+i)/60;
            i = p.edge_list(k,2);
            i_ = i_ + 1;
            if i == p.nn+1
                all_nodes_found = 1;
                P(i_,:) = [p.nn+1 0 0 pp(p.ne+p.ntc+2*p.nn+1)/60];
            end
        end
    end
end
zc = 1;
while zc < Lp && P(zc,2) == 0
    zc = zc + 1;
end
P = P(zc:end,:);

end