%% HW 5
%% Q1.b)
disp("Question 1.b)");
k = repmat(50, 4, 1);
C = [0 92 50 56
     92 0 80 74
     50 80 0 18
     56 74 18 0];
[y, TC] = ufladd(k, C);
fprintf(['Locations of the new facilities are at EF%d and EF%d,'...
    'and the total cost is %d.\n'], y(1), y(2), TC);
%% c)
disp("c)");
clear mp
mp = Milp('UFLADD');
[n m] = size(C);
kn = iff(isscalar(k),repmat(k,1,n),k(:)');
mp.addobj('min',kn,C)
for j = 1:m
   mp.addcstr(0,{':',j},'=',1);
end
for i = 1:n
   mp.addcstr({m,{i}},'>=',{i,':'});
end
mp.addub(1,1);
mp.addctype('B','C');

ilp = mp.milp2ilp;
[x,TC] = intlinprog(ilp{:});
x = mp.namesolution(x);
y = find(round(x.kn));
fprintf(['Locations of the new facilities are at EF%d, EF%d,'...
    'and EF%d, and the total cost is %d.\n'], y(1), y(2), y(3), int16(TC));
%% Q2
disp("Question 2)");
clear all
%% Input spreadsheet as table data
T1 = readtable('HW5data.xlsx','Sheet', 1);
T2 = readtable('HW5data.xlsx','Sheet', 2);
S1 = table2struct(T1);
S2 = table2struct(T2);
plts = uscity(mand({S2.City}, uscity('Name'),{S2.State}, uscity('ST')));
f = [S1.Demand]';
TPC = [S2.ProdCost]';
TDC = [S2.DistCost]';
%% Plot plants
makemap(plts.XY);
pplot(plts.XY,'rv');
pplot(plts.XY, plts.Name);
%% Get aggregate demand points
P = uszip5('XY', mor(T1.Zip,uszip5('Code5')));
pplot(P, 'g.')
%% Allocate plants demand from distributors based on distance
idx = zeros(1,length(f));
fDC = zeros(1,length(T2.City));
for i = 1:length(f)
   idx(i) = argmin(dists(P(i,:), plts.XY));
   fDC(idx(i)) = sum(fDC(idx(i)) + f(i));
end
% sum(fDC) - sum(T1.Demand) % checking if demand is correct

%% Nominal transport rate
length([T1.Zip]) == length(unique([T1.Zip]))
D = dists(plts.XY,P, 'mi');
F = sparse(argmin(D,1),1:length(S1),f);
r = sum(TDC)/sum(sum(F.*D));  
%% Est fixed cost
x = sum(F,2);
y = TPC;
yest = @(x,p) p(1) + p(2)*x;
fh = @(p) sum((y - yest(x,p)).^2);
ab = fminsearch(fh,[0 1]);
k = ab(1), cp = ab(2)
plot(x,y,'r.');
hold on
fplot(@(x) yest(x, ab), [0 max(x)],'k-'), hold off
%% Cost matrix for UFL
D = 1.2*dists([P; plts.XY],P, 'mi');
C = r*(f(:)'.*D);
%% UFL 
[y,TC,X] = ufl(k,C);
nNF = length(y)
TDC_new = sum(sum(C.*X))
TC_new = length(y)*k + TDC_new
TDC_orig = sum(TDC)
TC_orig = k * size(plts.XY,1) + TDC_orig
100*TDC_new/TDC_orig
100*TC_new/TC_orig
%% Plot new DCs
NFsites = [P; plts.XY];
pplot(NFsites(y,:),'bo');
i = int2str([1:length(y)]');
pplot(NFsites(y,:), cellstr(i));
%% Decision to make
clear ind 
clear i x
for i=1:nNF
    x{i} = struct('NF', {i},...
        'Distance', {min(1.2*dists(NFsites(y(i),:),plts.XY,'mi'))},...
        'ClosestPlant',{plts.Name(argmin(1.2*dists(NFsites(y(i),:),...
        plts.XY,'mi')))}); 
    if x{i}.Distance > 100
        fprintf('NF at the site %d should be constructed.\n', x{i}.NF)
    end
end
for i=1:length(plts.Name)
    if min(1.2*dists(plts.XY(i,:), NFsites(y,:), 'mi')) > 100
        ind(i) = i;
        fprintf('EF in %s should be closed.\n', char(plts.Name(i)))
    end
end
pplot(plts.XY(nonzeros(ind'),:),'kx','MarkerSize',12)
%% Q3
disp("Question 3)");
clear all
close all
s = uscity(strcmp('NC',uscity('ST')) & uscity('Pop') > 20000);
rmax = 30;
D = dists(s.XY,s.XY,'mi');
D = D + sqrt(s.LandArea/pi);
c = ones(1,size(D,2));
A = false(size(D));
A(D <= rmax) = true;
is0 = ~any(A,2);
A(is0,:) = [];

mp = Milp('Set Covering');
mp.addobj('min',c);
mp.addcstr(A,'>=',1);
mp.addctype('B');

ilp = mp.milp2ilp;
x = intlinprog(ilp{:});
nNF = sum(x);

idx = find(x);
makemap(s.XY)
pplot(s.XY,'r.')
pplot(s.XY(idx,:),'go')
pplot(s.XY(idx,:),s.Name(idx))
fprintf('New transmitters should be installed at %s.\n',...
    join(string(s.Name(idx)), ', '))
PopCov = 100*sum(s.Pop)/sum(nccity('Pop'));
fprintf('These transmitters will cover %.2f%% of NC population.\n',PopCov)
%% Q4
disp("Question 4)");
clear all
clear mp
x = uszip5(mor({'NC'},uszip5('ST')) & uszip5('Pop') > 30000)
P = x.XY;
a = x.LandArea;
v = x.Pop';
dafh = @(XY1,a1,XY2,a2) max(1.2*dists(XY1,XY2,'mi'),...
   0.675*max(sqrt(a1(:)),sqrt(a2(:)')));
Da = dafh(P,a,P,a);  % Area-adjusted distances
f = x.Pop';                 % person
r = 1/10000;                % $/person-mi
C = r*Da.*f;                % $
k = 500;  % repmat(500,1,size(C,1));
TCold = 1.2102e+04;
%% MILP model
clear mp
mp = Milp('UFL')
mp.Model
[n m] = size(C)
M = 1:m;
V = 400000;
kn = iff(isscalar(k),repmat(k,1,n),k(:)');  % expand if k is constant value
mp.addobj('min',kn,C)  % min sum_i(ki*yi) + sum_i(sum_j(cij*xij))
for j = 1:m
   mp.addcstr(0,{':',j},'=',1)   % sum_i(xij) = 1
end
for i = 1:n
   mp.addcstr({m,{i}},'>=',{i,':'})  % m*yi >= sum_j(xij)  (weak form.)
end
for i = M
   mp.addcstr({V,{i}},'>=',{v,{i,':'}})   % V*yi >= sum_j(vj*xij)
end
mp.addub(1,1)
mp.addctype('B','C')        
mp.Model
%% Solution
ilp = mp.milp2ilp;
[x,TC,exitflag,output] = intlinprog(ilp{:});
x = mp.namesolution(x) 
nNF = sum(x.kn);
TCchng = 100*(TC-TCold)/TCold;
fprintf('Number of facilities needed is %d.\n',nNF)
fprintf('Change in total cost is +%.2f%%.\n',TCchng)
