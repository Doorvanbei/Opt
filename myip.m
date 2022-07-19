function [x,v,f] = myip(c,iind,Al,bl,A,b)
% MILP solver
% doorvanbei
% 20220711
pvn = length(c); % primal variable number
cutnum = 0; % count how many gomory cuts added in the process
ub = inf; % ub for IP opt val
T = 1; % T is nonempty at first, root node is node 1
Tind = []; % record tree info
cutpo = {}; % global space storing gomory cuts
% cutpo{r,1} = All(Tree_r), cutpo{r,2} = bll(Tree_r), r = 1,2,...
while ~isempty(T) % tree has not been fully searched
    currNode = T(1); % for the ease of Debugging
    T(T==currNode) = []; % firstly delete this node (then deal with this LP completely)
    All = []; % store local GomoryMI cuts at currNode
    bll = []; % 'll': lessEqual, local
    if currNode > 1 % restrict somewhere
        % decode branching tree: which variables to be branched and how to
        xbind = Tind(currNode,1);
        xbval = Tind(currNode,2);
        xbdir = Tind(currNode,3);
        lastNode = Tind(currNode,4);
        while lastNode~=1 % not till the root tree
            xbind = [xbind Tind(lastNode,1)];
            xbval = [xbval Tind(lastNode,2)];
            xbdir = [xbdir Tind(lastNode,3)];
            lastNode = Tind(lastNode,4);
        end
        % branching constraints are all variable bound constraints
        l = length(xbind);
        xlb = zeros(pvn,1); % initialize
        xub = inf*ones(pvn,1); % initialize
        for i = 1:l
            if xbdir(i) % do the right branch
                nlb = ceil(xbval(i));
                if nlb > xlb(xbind(i)) % new lower bound is better than current
                    xlb(xbind(i)) = nlb; % update it
                end
            else % do the left branch
                nub = floor(xbval(i));
                if nub < xub(xbind(i))
                    xub(xbind(i)) = nub;
                end
            end
        end
        eqi = xlb == xub; 
        neqi = ~eqi; % after solving reduced sys x_red, write x(eqi) = x_fixed, x(neqi) = x_red;
        boundeqnum = sum(eqi); % how many 'xi=z0'
        % write (the rest) variable bound constraints to Al * x <= bl
        ltmp = find( xlb>0 & neqi );
        utmp = find( xub<inf & neqi );
        ltl = length(ltmp);
        utl = length(utmp);
        blb = zeros(ltl+utl,1); % store ineq constraints corr to branching constraints
        Alb = zeros(ltl+utl,pvn);
        for i = 1:utl
            Alb(i,utmp(i)) = 1;
            blb(i) = xub(utmp(i));
        end
        for i = 1:ltl
            Alb(utl+i,ltmp(i)) = -1;
            blb(utl+i) = -xlb(ltmp(i));
        end
        % col reduction
        x_fixed = xlb(eqi);
        if ~isempty(b)
            local_b = b - A(:,eqi) * x_fixed;
            local_A = A(:,neqi);
        else
            local_b = [];
            local_A = [];
        end
        if ~isempty(bl)
            local_bl = bl - Al(:,eqi) * x_fixed;
            local_Al = Al(:,neqi);
        else
            local_bl = [];
            local_Al = [];
        end
        local_c = c(:,neqi);
        if ~isempty(Alb)
            Alb = Alb(:,neqi);
        end
        earlyflag = 0;
        while 1 % make formulation of this node better, try to set earlyflag = 1
            [x_red,~,f,info] = glp(local_c,[local_Al;Alb;All],[local_bl;blb;bll],local_A,local_b);
            if f == -2 % currNode is infeas: actually this is not a leave node
                earlyflag = 1;
                break
            elseif f == -3
                error('In mypureipbac: Unexpected Err!')
            end
            x_red = x_red(1:(pvn-boundeqnum)); % LP relax is normal
            x = zeros(pvn,1);
            x(eqi) = x_fixed;
            x(neqi) = x_red;
            v = c * x;
            if v >= ub % pruned by bound
                earlyflag = 1;
                break
            end
            if all(isint(x(iind))) % int attained at curr Node
                ub = v;
                xOpt = x; % current opt S-feas point.
                earlyflag = 1;
                break
            end
            % try to add cuts before directly go to branching
            [lhs,rhs,failflag] = gencut([local_Al;Alb;All],[local_bl;blb;bll],info{3},info{2},info{1},pvn-boundeqnum,iind);
            if failflag % cut gened is not proper
                cutpo{currNode,1} = All;
                cutpo{currNode,2} = bll;
                break % stop adding more cuts and then go to branching
            else % store this cut in ax<=b form
                All = [All;lhs];
                bll = [bll;rhs];
                cutnum = cutnum + 1;
            end
        end
        if earlyflag
            continue
        end
        % do tree branching
        [~,bi] = max(abs(x(iind) - round(x(iind))));
        bi = iind(bi);
        % 1st row of Tind ==> tree (LP problem) number = 1 (1 means the root)
        % 1st col of Tind: bi = 2 ==> do branching to x_2
        % 2nd : save the fractional val to be branched
        % 3rd : 0 - do the left branch, 1 - do the right branch
        % 4th : 0 ==> the father tree of 'the tree corresponding to this row' is 0
        Tind = [Tind;bi,x(bi),0,currNode];
        Tind = [Tind;bi,x(bi),1,currNode];
        [mT,~] = size(Tind);
        T = [T mT-1 mT];
    else % LP0 relax
        % Gomory cuts gen.
        while 1 % We gen suitable All,bll in this loop
            [x,v,f,info] = glp(c,[Al;All],[bl;bll],A,b);
            if f ~= 1
                disp('In mypureipbac: root problem not normal, see exit f below.')
                disp(f)
                return
            end
            x = x(1:pvn); % do not consider injecting variables.
            if all(isint(x(iind))) % Optimality attained at root node.
                v = c * x;
                f = 1;
                return
            end
            % add some gomory cuts.
            [lhs,rhs,failflag] = gencut([Al;All],[bl;bll],info{3},info{2},info{1},pvn,iind); % pvn = length(c_local)
            if failflag % cut gened is not proper
                cutpo{currNode,1} = All;
                cutpo{currNode,2} = bll;
                break
            else % store this cut in ax<=b form
                All = [All;lhs];
                bll = [bll;rhs];
                cutnum = cutnum + 1;
            end
        end
        % tree branching
        [~,bi] = max(abs(x(iind) - round(x(iind))));
        bi = iind(bi);
        Tind = [0,0,0,0]; % no info for the root Tree
        Tind = [Tind;bi,x(bi),0,currNode]; % create Tree 2
        Tind = [Tind;bi,x(bi),1,currNode]; % create Tree 3
        [mT,~] = size(Tind); % mT = 2
        T = [T mT-1 mT]; % T = [2,3]
    end
end
% disp('cutnum:')
% disp(cutnum)
x = xOpt;
v = c * xOpt;
f = 1;

