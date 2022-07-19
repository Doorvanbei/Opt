function [x,v,f,info] = glp(c,Al,bl,A,b)
% please use this function only when A,Aeq are both nonempty.
% mix cx | Al x <= bl, Ax = b, x>=0, where all params are INT excluding c.
    [m,~] = size(Al);
    [n,~] = size(A);
    [x,v,f,info] = mylp([c zeros(1,m)],[A zeros(n,m);Al eye(m)],[b;bl]);
end