function [category,theta,nc,ntc] =...
    get_categories(nn,aa,rr,user_distributions)
La = length(aa);
Lr = length(rr);
C = zeros(La,Lr,nn); % categories (full version)
nu = length(user_distributions(:,1));
% idx(aik) | idx(rik) | i
for ja = 1:La
    for jr = 1:Lr
        for k = 1:nu
            i = user_distributions(k,2);
            ak = user_distributions(k,3);
            rk = user_distributions(k,4);
            if ak <= aa(ja) && aa(ja) < ak+rk && rr(jr) <= ak+rk-aa(ja)
                C(ja,jr,i) = C(ja,jr,i) + 1;
            end
        end
    end
end
nc = zeros(nn,1); % number of categories for each node
%fprintf('Categories (idx(aik) | idx(rik) | i):\n\n')
category = cell(nn,1); % categories (sparse version)
% (idx(ail), idx(ril)) -> nil
theta = cell(nn,nu);
for i = 1:nn
    [I,J,K] = find(C(:,:,i));
    nc(i) = length(K);
    for k = 1:nu
        if i == user_distributions(k,2)
            theta{i,k} = zeros(nc(i),1);
        end
    end
    for l = 1:nc(i)
        I(l) = aa(I(l));
        J(l) = rr(J(l));
        for k = 1:nu
            if i == user_distributions(k,2)
                ak = user_distributions(k,3);
                rk = user_distributions(k,4);
                if ak <= I(l) && I(l) < ak+rk && J(l) <= ak+rk-I(l)
                    theta{i,k}(l) = 1;
                end
            end
        end
    end
    category{i} = [I J K];
    %fprintf(['Node ' num2str(i) '\n'])
    %display(category{i})
end
ntc = sum(nc);

end
