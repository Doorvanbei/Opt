function [lhs,rhs,failflag] = gencut(Al,bl,Ta,dind,bind,pvn,iind)
% Al: initial formulation of Al
% bl: initial formulation of bl
% pvn = 3; % primal variable number, excluding those injecting variables by
%       ineq constraints.
% dind = [2 5]; % => a cut is of form:  c_ * x2 + c_ * x5 >= 1

temp = Ta(1,2:end); % a template for coeff of all variables (inj included)
temp(:) = 0;
tmpi = temp;
tmpb = temp;
tmpb(bind) = 1;
tmpi(iind) = 1;
biind = find(tmpb & tmpi);
l = length(biind);
rnumber = biind;
for i = 1:l
    rnumber(i) = find(bind==biind(i));
end
final_b = Ta(rnumber,end);
[~,r] = max(abs(final_b - round(final_b))); % use row r
rhs = final_b(r);
r = rnumber(r); % gen cut from row r
f0 = rhs - floor(rhs); % 0<f0<1

NIind = find(~tmpb &  tmpi); % []
NCind = find(~tmpb & ~tmpi); % [4,5,6]

% dind is equal to ind1 | ind2 | ind3 | ind4
f = Ta(r,NIind);
f = f - floor(f);
ind1 = NIind(f<=f0);
ind2 = NIind(f>f0);

a = Ta(r,NCind);
ind3 = NCind(a>=0);
ind4 = NCind(a<0);

nr = temp;
nr(ind3) = Ta(r,ind3)/f0;
nr(ind4) = -Ta(r,ind4)/(1-f0);
f = Ta(r,ind1);
f = f - floor(f);
nr(ind1) = f/f0;
f = Ta(r,ind2);
f = f - floor(f);
nr(ind2) = (1-f)/(1-f0);
nr = nr(dind);

if isempty(bl)
    lhs = zeros(1,pvn);
    lhs(dind) = -nr;
    rhs = -1;
else
    r = zeros(1,pvn+length(bl));
    r(dind) = nr; % r * [pv;injv] >= 1
    rL = r(1:pvn);
    rR = r((1+pvn):end); % rL * pv >= 1 - rR * injv
    % injv = bl - Al * pv
    lhs = rR * Al - rL;
    rhs = rR * bl - 1;
end
failflag = 0;
if ~iscutok(lhs,rhs,25)
    failflag = 1;
end

