% the main .m file used to test prisimp.m
% doorvanbei
% 20220702

while 1
epsi = 1e-6;
cols = 5;
rows = 8;
A0 = (randi(3,[rows,cols]) ./ (randi(5,[rows,cols])-2.3)) .* (randi(2,[rows,cols])-1);
b0 = randi(6,[rows,1])-3;
c0 = randi(6,[1,cols])-3;

% in this case：matlab solver say 'unbound', my solver say 'infeasible'
% using fourier elim --> I conclude 'infeasible'
% A0 =[   -1.5385         0    0.3704   -3.3333   -0.7692
%          0   -0.7692         0         0         0
%     0.7407    4.2857         0         0         0
%     0.5882    0.5882         0         0    0.7407
%    -2.3077    1.1765         0         0         0
%          0         0         0   -0.7692    1.1765
%     0.7407         0         0         0         0
%     0.5882    1.7647         0         0   -3.3333];
% 
% b0 =[    -2
%     -1
%      2
%      1
%     -1
%      2
%      2
%     -2];
% 
% c0 =[     1     2    -1    -1     3];

% this case should be unbound but matlab solver says 'opt attained.'
% A0 = [         0         0   -1.5385         0         0
%     2.8571   -1.5385   -1.5385         0    0.7407
%    -2.3077   -3.3333         0         0         0
%     2.8571   -2.3077    1.7647         0         0
%    -0.7692   -2.3077         0         0         0
%    -2.3077         0   -6.6667         0         0
%    -3.3333    0.5882   -3.3333    1.1765   -0.7692
%    -0.7692         0   -2.3077   -1.5385         0];
% 
% b0 =[    -1
%     -1
%      0
%      2
%     -1
%     -1
%      3
%      1];
% 
% c0 = [     0     3    -2    -2     2];
% 
% x1 = [    0.1978
%     0.3674
%     0.6500
%     4.7683
%          0];

if rank(A0)<cols
    continue
end

A = A0;
b = b0;
c = c0;

[x0,v0,f0] = linprog(c,A,b,[],[],[0;0;0;0;0],inf*[1;1;1;1;1]); % use matlab solver

[m,n] = size(A); 
c = [c zeros(1,m)];
A = [A eye(m)];
[m,n] = size(A); 
if rank(A) < m
    continue
end
[x,v,f] = prisimp(c,A,b); % use my solver

if f0 == -2
    if f ~= f0
        error('infeasible err!')
    end
    continue
end

if f0 == -3
    if f == -2
        continue % omit this case temporarily
    elseif f~=f0
        error('unbound err!')
    end
    continue
end

if f0 == 1
    if f~=f0
        error('normal: type err!')
    else
%         v
%         v0
%         fprintf('val distance %f\n',abs(v-v0));
        if abs(v-v0) > epsi
            error('normal: opt val v has a dist with std v0!')
        end
        p = A * x - b;
        if any(p>epsi)
            error('normal: prisimp.m return an infeasible opt solution!')
        end
    end
end

end
 
% A = [3 1 -2 1 -3;1 3 0 -1 -1];
% b = [6;2];
% c = [-19 12 2 6 -19];
% 
