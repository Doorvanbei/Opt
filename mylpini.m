function [bind,info] = mylpini(c,A,b)
    % the artificial LP for initializing primal simplex method
    % doorvanbei
    % 20220704
    epsi = 1e-6;
    info = [];
    [m,n0] = size(A);
    % starting: make sure rhs vector b >= 0
    tmp = b<0;
    A(tmp,:) = -A(tmp,:);
    b(tmp) = -b(tmp);
    % add artificial variables
    A = [A eye(m)];
    c(:) = 0;
    c = [c ones(1,m)];
    [~,n] = size(A);
    nind = 1:n;
    bind = nind((end-(m-1)):end); % initial Basis
    % simplex method for solving sys min cx | Ax=b,x>=0
    while 1
        dind = nind;
        dind(bind) = [];
        B = A(:,bind);
        a0 = B\b; % >=0 : feasible initial basis, >0 : nondegen
        if any(a0<-epsi) % make sure that primal side is always feasible
            error('Initial Basis Not Feasible!')
        end
        cb = c(bind);
%         val = cb * a0; % intermediate value, for debugging
        cd = c(dind);
        y = cb/B; % shadow price
        D = A(:,dind);
        rd = cd - y * D;
        if ~any(rd<-epsi) % Once dual side become feasible, adding that primal val = dual val, opt attained.
            if cb * a0 > epsi % if opt val > 0
                x = zeros(n,1);
                x(bind) = a0; % save this primal opt solu
                info = {};
                info{1} = c;
                info{2} = A;
                info{3} = b; % save Artificial LP sys
                info{4} = x; % primal opt solu. check1: is it primal feas? (A*x-b=0, x>=0)
                info{5} = y; % dual opt solu. check2: is y dual feas? (y*A-c <= 0)
                % check3: Does pval = dval?  c * x = y * b (should > 0)
                bind = []; % region Ax=b, x>=0 is infeasible
            elseif any(bind>n0)
               tmp = bind>n0;               
               currBind = bind;
               currBind(tmp) = []; % [1 4]
               restBind = 1:n0;
               restBind(currBind) = []; % [2 3 5]
               Acurr = A(:,currBind);
               cho = nchoosek(restBind,sum(tmp));
               [r,~] = size(cho);
               f = 1;
               for i = 1:r
                   Arest = A(:,cho(i,:));
                   if rank([Acurr Arest]) == m
                       bind = [currBind cho(i,:)];
                       f = 0;
                       break
                   end
               end
               if f
                   error('mylpini fail to give a initial BFS!');
               end

            end
            break
        else % curr BFS can be strictly improved
            [~,i] = min(rd); % entering basis vector in col dind(i)
            ae = B\D(:,i);
            r = a0./ae;
            r(ae<epsi) = inf; % degeneracy included
            [~,o] = min(r);
            bind(o) = dind(i);
        end
    end
end
