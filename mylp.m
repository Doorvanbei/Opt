function [x,v,f,info] = mylp(c,A,b)
	% min cx | Ax = b, x >= 0
    % return optSolu,optVal,exitFlag: 1 for normal,-2 for infeas,-3 for unb
    % algorithm: primal simplex method
    % doorvanbei
    % 20220704
    epsi = 1e-6;
    %info = [];
    [m,n] = size(A);
    if ~(m<n && rank(A) == m && all(ones(1,m)*abs(A)>0))
        error('Input LP param NOT in std form! Make sure that A is rectangle,full-row-rank and contains no 0-vec columns')
    end
    nind = 1:n;
    [bind,info] = mylpini(c,A,b); % initialization using artificial LP
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
%         val = cb * a0; % intermediate value, for debugging
        cd = c(dind);
        y = cb/B; % shadow price
        D = A(:,dind);
        rd = cd - y * D;
        if ~any(rd<-epsi) % opt cond: 1,p side feas 2,d side feas 3,val eq at 2 sides.
            T = [B\[A b]];
            ct = [c 0];
            ct = ct - cb * T;
            info = {};
            info{1} = bind; % opt basis index
            info{2} = dind;
            info{3} = [T;ct]; % opt Tableau
            v = cb * a0; % also eq2 'y*b'
            x = zeros(n,1);
            x(bind) = a0;
            f = 1;
            return
        else % curr BFS can be strictly improved
            [~,i] = min(rd); % entering basis vector in col dind(i)
            ae = B\D(:,i);
            if ~any(ae>epsi) % unbound
                dind(i) = [];
                le = length(dind);
                Aray2 = zeros(le,n);
                for i = 1:le
                    Aray2(i,dind(i)) = 1;
                end
                Aray3 = zeros(1,n);
                Aray3(1) = 1;
                Aeq = [A;Aray2;Aray3];
                cnt = 0;
                ra = rank(Aeq);
                while ra<n
                    Aeq(n,:) = circshift(Aeq(n,:),1); % 1not, 2not, try n
                    cnt = cnt + 1;
                    if cnt >= n
                        error('In case f=-3: no valid ray exists!')
                    end
                    ra = rank(Aeq);
                end
                info = Aeq\[zeros(m+le,1);1]; % output extRay corr. to final vertex.
                f = -3; % unbounded flag
                x = zeros(n,1);
                x(bind) = a0; % final vertex
                v = cb * a0; % final feas val
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


















