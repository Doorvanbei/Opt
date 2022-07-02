function bind = innerprimsimp(c,A,b)
    % the artificial LP for initializing primal simplex method
    % doorvanbei
    % 20220702
    epsi = 1e-6;
    [m,~] = size(A);
    % starting: make sure rhs vector b >= 0
    tmp = b<0;
    A(tmp,:) = -A(tmp,:);
    b(tmp) = -b(tmp);
    % add artificial variables
    A = [A eye(m)];
    c(:) = 0;
    c = [c ones(1,m)];
    [m,n] = size(A);
    nind = 1:n;
    bind = nind((end-(m-1)):end); % initial Basis
    while 1
        dind = nind;
        dind(bind) = [];
        B = A(:,bind);
        a0 = B\b; % >=0 : feasible initial basis, >0 : nondegen
        if any(a0<-epsi) % make sure that primal side is always feasible
            error('Initial Basis Not Feasible!')
        end
        cb = c(bind);
%         val = cb * a0 % intermediate value, for debugging
        cd = c(dind);
        y = cb/B; % shadow price
        D = A(:,dind);
        rd = cd - y * D;
        if ~any(rd<-epsi) % Once dual side become feasible, adding that primal val = dual val, opt attained.
            if cb * a0 > 0 % if opt val > 0
                bind = []; % region Ax=b, x>=0 is infeasible.
            end
            return
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








%             rdred = rd; % save original rd
%             [minuv,i] = min(rdred); % entering basis vector in col dind(i)
%             while minuv < 0
%                 ae = B\D(:,i)
%                 r = a0./ae
%                 r(ae<epsi) = inf
%                 [~,o] = min(r);
%                 bindtmp = bind;
%                 bindtmp(o) = dind(i);
%                 qB = A(:,bindtmp);
%                 if rank(qB) < m
%                     disp('(((((((((((((((((((((((((((((((((())))))')
%                     rdred(i) = inf; % modify rdred
%                     [minuv,i] = min(rdred); % gen a new quasi-enter index i
%                     continue
%                 else % regular case
%                     bind(o) = dind(i);
%                     break
%                 end
%             end
%             if minuv >= 0
%                 error('curr Basis not opt, but cannot find entering vec!')
%             end