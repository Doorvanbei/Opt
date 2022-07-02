function [x,v,f] = prisimp(c,A,b)
	% min cx | Ax = b, x >= 0
    % return optSolu,optVal,exitFlag: 1 for normal,-2 for infeas,-3 for unb
    % algorithm: primal simplex method
    % doorvanbei
    % 20220702
    epsi = 1e-6;
    [m,n] = size(A);
    if ~(m<n && rank(A) == m)
        error('Input LP param NOT in std form!')
    end
    nind = 1:n;
    bind = innerprimsimp(c,A,b); % initialization using artificial LP
    if isempty(bind)
        f = -2; % LP infeasible
        x = [];
        v = inf;
        return
    end
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
        if ~any(rd<-epsi) % opt cond: 1,p side feas 2,d side feas 3,val eq at 2 sides.
            v = cb * a0; % also eq2 'y*b'
            x = zeros(n,1);
            x(bind) = a0;
            f = 1;
            return
        else % curr BFS can be strictly improved
            [~,i] = min(rd); % entering basis vector in col dind(i)
            ae = B\D(:,i);
            if ~any(ae>epsi) % unbound
                f = -3;
                x = inf;
                v = -inf;
                return
            else
                r = a0./ae;
                r(ae<epsi) = inf;
                [~,o] = min(r);
                bind(o) = dind(i);
            end
        end
    end
end


















