%% HW7
%% Q1 input
clear
disp('Question 1')
s = [20548 26149 36317];
dem = [40 55 35 70 25];
c = [30669 38339 30732 23830 23154];
Sup = uszip5(mand(s, uszip5('Code5')));
Cust = uszip5(mand(c, uszip5('Code5')));
%% a)
disp('a)')
D = dists(Sup.XY, Cust.XY, 'mi');
% mdisp(D)
% argmin(D, 1)
F = full(sparse(argmin(D,1),1:length(dem),dem))
TCa = D.*F;
TCa = sum(sum(D.*F));
for i=1:length(s)
    if any(F(i,:))
        fprintf(['DC%d should supply %s products to customers in %s, '...
        'respectively\n'], i, int2str(nonzeros(F(i,:))'),...
        int2str(c(find(any(F(i,:),1)))));
    end
end
fprintf('Total ton-miles is %.2f.\n', TCa)
%% b)
disp('b)')
sup = [60 90 80]';
% mdisp([D sup; dem 0])
[F,TCb] = trans(D,sup,dem);
TCb
fprintf('The change in total ton-miles is %.2f.\n', TCb-TCa);
%% c) 
disp('c)')
%% Create MCNF inputs
IJC = lev2list(D);
IJCU = [IJC repmat(30, length(IJC),1)];
s = [sup' -dem];
%% solve
lp = mcnf2lp(IJCU,s);
[x,TCc,XFlg,out] = lplog(lp{:});
[f,TCc,nf] = lp2mcnf(x,IJC,s);
%% report
IJF = [IJC(:,[1 2]) f];
AF = list2adj(IJF);
F = adj2lev(AF,size(D))
TCc
fprintf('The change in total ton-miles is %.2f.\n', TCc-TCb);
%% d)
disp('d)')
%% input
p = [28124 27325 37421 27513];
Plts = uszip5(mand(p, uszip5('Code5')));
%% Create MCNF inputs
C23 = D;
C12 = dists(Plts.XY, Sup.XY, 'mi');
W = lev2adj(C12,C23);
IJC = adj2list(W);
s = [repmat(60, 4,1)'  zeros(3,1)' -dem];
    
%% solve
lp2 = mcnf2lp(IJC,s);
[x,TCd] = lplog(lp2{:});
[f,TC,nf] = lp2mcnf(x,IJC,s);
%% Report
IJF = [IJC(:,[1 2]) f];
AF = list2adj(IJF);
[F12,F23] = adj2lev(AF,[4 3 5])

for i=1:4
    fprintf(['Plant%d should supply %s products to DCs in %s, '...
    'respectively\n'], i, int2str(nonzeros(F12(i,:))'),...
    int2str(p(find(any(F12(i,:),1)))));
end
for i=1:3
    if any(F23(i,:))
        fprintf(['DC%d should supply %s products to customers in %s, '...
        'respectively\n'], i, int2str(nonzeros(F23(i,:))'),...
        int2str(c(find(any(F23(i,:),1)))));
    end
end

%% Q2 input
clear
disp('Question 2b)')
%% Create data 
IJD = [
    1 -2 14
    1 -6 17
    1 -8 13
    1 -9 1
    1 -10 16
    2 -3 2
    2 -4 6
    2 -7 16
    2 -9 10
    2 -10 7
    3 -4 5
    3 -8 8
    3 -9 9
    4 -5 14
    4 -7 14
    5 -6 7
    5 -7 1
    5 -10 16
    6 -8 19
    6 -10 10
    7 -10 1
    8 -9 10];
%% Dijkstra's algorithm
[d,p] = dijkdemo(list2adj(IJD),3,6)
%% Q3 input
clear
disp('Question 3')
Ral = [-78.701389 35.7725];
Atl = [-84.39 33.771944];
XY1 = [Ral; Atl];
%% Get road network
expansionAroundXY = 0.1;
[XY2,IJD,isXY,isIJD] = subgraph(usrdnode('XY'),...
   isinrect(usrdnode('XY'),boundrect(XY1,expansionAroundXY)),...
   usrdlink('IJD'));
%% Label type of road
s = usrdlink(isIJD);
isI = s.Type == 'I';         % Interstate highways
isIR = isI & s.Urban == ' '; % Rural Interstate highways
isIU = isI & ~isIR;          % Urban Interstate highways
isR = s.Urban == ' ' & ~isI; % Rural non-Interstate roads
isU = ~isI & ~isR;           % Urban non-Interstate roads
%% Plot roads
makemap(XY2,0.03)  % 3% expansion
h = [];  % Keep handle to each plot for legend
h = [h pplot(IJD(isR,:),XY2,'r-','DisplayName','Rural Roads')];
h = [h pplot(IJD(isU,:),XY2,'k-','DisplayName','Urban Roads')];
h = [h pplot(IJD(isI,:),XY2,'c-','DisplayName','Interstate Roads')];
%% Add connector roads from cities to road network
[IJD11,IJD12,IJD22] = addconnector(XY1,XY2,IJD);
h = [h pplot(IJD12,[XY1; XY2],'b-','DisplayName','Connector Roads')];
h = [h pplot(XY1,'go','DisplayName','Destinations')];
%% Convert road distances to travel times (needs to be after ADDCONNECTOR)
v.IR = 75;  % Rural Interstate highways average speed (mph)
v.IU = 65;  % Urban Interstate highways average speed (mph)
v.R = 50;   % Rural non-Interstate roads average speed (mph)
v.U = 25;   % Urban non-Interstate roads average speed (mph)
v.C = 20;   % Facility to road connector average speed (mph)

IJT = IJD;
IJT(isIR,3) = IJD(isIR,3)/v.IR;
IJT(isIU,3) = IJD(isIU,3)/v.IU;
IJT(isR,3) = IJD(isR,3)/v.R;
IJT(isU,3) = IJD(isU,3)/v.U;

IJT22 = IJD22;                % road to road
IJT22(:,3) = IJT(:,3);
IJT12 = IJD12;                % facility to road
IJT12(:,3) = IJD12(:,3)/v.C;  % (IJD11 facility to facility arcs ignored)

%% Find shortest path
[~,P] = dijk(list2adj([IJT12; IJT22]),1:2);
[T, p] = dijk(list2adj([IJT12; IJT22]),1,2);

%% Distance of shortest time route

W = list2adj([IJD12; IJD22]);
D = locTC(pred2path(P,1,2),W);

h = [h ...
   pplot({p},[XY1;XY2],'y-','LineWidth',2,'DisplayName','Shortest Path')];
title(sprintf(['From home in Raleigh to colleague in GaTech:\n'...
    'Distance %.2f mi, Time = %d hr %d min'],D,floor(T),round(60*(T-floor(T)))));
legend(h),shg
fprintf(['From home in Raleigh to colleague in GaTech: Distance %.2f mi,'...
    'Time = %d hr %d min.\n'],D,floor(T),round(60*(T-floor(T))));
%% Q4 input
clear
disp('Question 4')
T = 26;
rng(1964);
D = round([gamrnd(6,4,T,1) gamrnd(4,3,T,1)]);

K = [60 50;
     55 45;
     50 35];

Cp = [12  20;     
      75 130;
      35  60];
h = 0.4/365.25*7;
Ci = cumsum(Cp,1)*h      % inventory cost of product g for stage m ($/ton)

Cs = [400 600;         % stage-m product-g setup cost ($)
       90  110;
       50  60];
yinit = [0  0;           % initial product g inventory at stage m (ton)
         0  0;
         0 sum(D(1:2,2))];
yfinal = zeros(3,2);     % final product g inventory at stage m (ton)
k0 = [1 0;               % initial setup at stage m for product g
      1 0;
      1 0];
M = size(K,1);           % number of production stages = 3
T = size(D,1);           % number of periods of production = 6
G = size(K,2);           % number of products produced = 2
%% Create MILP model
Cp = reshape(repmat(Cp,[T 1 1]),M,T,G);     % create M x T x G array (3-D)
Ci = reshape(repmat(Ci,[T+1 1 1]),M,T+1,G); % create M x (T+1) x G array
Ci(:,1,:) = 0;   % intital inventory cost already accounted for last period
Cs = reshape(repmat(Cs,[T 1 1]),M,T,G);     % create M x T x G array
mp = Milp('PPlan');
mp.addobj('min',Cp,Ci,Cs,zeros(M,T,G));     % zeros(M,T,G) dummy array for k
for g = 1:G
   for t = 1:T
      for m = 1:M-1
         mp.addcstr({[1 -1],{[m m+1],t,g}},{[1 -1],{m,[t t+1],g}},0,0,'=',0)
      end
      mp.addcstr({M,t,g},{[1 -1],{M,[t t+1],g}},0,0,'=',D(t,g))
      for m = 1:M
         mp.addcstr({m,t,g},0,0,'<=',{K(m,g),{m,t,g}})
      end
   end
   for m = 1:M
      mp.addcstr(0,0,{-1,{m,1,g}},{m,1,g},'<=',k0(m,g))
      for t = 2:T
         mp.addcstr(0,0,{-1,{m,t,g}},{[1 -1],{m,[t t-1],g}},'<=',0)
      end
   end
end
for m = 1:M, for t = 1:T, mp.addcstr(0,0,0,{m,t,':'},'=',1), end, end
mp.addlb(0,horzcat(reshape(yinit,M,1,G),zeros(M,T-1,G),reshape(yfinal,M,1,G)),0,0);
mp.addub(Inf,horzcat(reshape(yinit,M,1,G),inf(M,T-1,G),reshape(yfinal,M,1,G)),1,1);
mp.addctype('C','C','B','B');

%% Solve using Gurobi
clear params
model = mp.milp2gb;
params.outputflag = 1;
result = gurobi(model, params);
x = mp.namesolution(result.x);
TC = result.objval
%% Report results

Fp = x.Cp;
Fi = x.Ci;
for g = 1:G
   mdisp(D(:,g)',[],[],['D' num2str(g)])
   mdisp(Fp(:,:,g),[],[],['Fp' num2str(g)])
   mdisp(Fi(:,:,g),[],[],['Fi' num2str(g)])
end
%% Q5 input
clear
disp('Question 5')
Q = [20;30];    % capacity of stage m in period t (ton)
D = [25 15 10 50 25 15]';  % demand in period t (ton)
Cp = [200; 800];           % production cost in stage m ($/ton)
Cf = [3000; 9000];          % fixed cost for stage m ($/ton)
h_ob = 0.08/(3/12);
h = (0.11+h_ob)/12;
Ci = cumsum(Cp+Cf./Q)*h;          % inventory cost for stage m ($/ton)

yinit = [5; 0];            % initial inventory of stage m (ton)
yfinal = [7; 4];            % final inventory of stage m (ton)
ymax = 30;
M = size(Q,1);             % number of production stages = 2
T = size(D,1);             % number of periods of production = 6
%% Create MILP model
Cp = repmat(Cp,1,T);     % create M x T array
Ci = repmat(Ci,1, T+1); % create M x (T+1) array
Ci(:,1,:) = 0;   % intital inventory cost already accounted
Cf = repmat(Cf,1,T);     % create M x T array

mp = Milp('PPlan');
mp.addobj('min',Cp,Ci,Cf);

for t = 1:T
   for m = 1:M-1
      mp.addcstr({[1 -1],{[m m+1],t}},{[1 -1],{m,[t t+1]}},0,'=',0)
   end
   mp.addcstr({M,t},{[1 -1],{M,[t t+1]}},0,'=',D(t))
   for m = 1:M
      mp.addcstr({m,t},0,'<=',{Q(m),{m,t}})
   end
end

mp.addlb(0,[yinit zeros(M,T-1) yfinal], 0)
mp.addub(Inf,[yinit repmat(ymax,M,T-1) yfinal], 2)  
mp.addctype('C','C','I')
%% Solve using Gurobi
ilp = mp.milp2ilp;
[x,TC,exitflag,output] = intlinprog(ilp{:});
TC, output
x = mp.namesolution(x);
%% Report results
Fp = x.Cp;
Fi = x.Ci;
Ff = x.Cf;
disp('Production plan is described below')
mdisp(D)
mdisp(Fp)
mdisp(Fi)
mdisp(Ff)
fprintf('The total cost with this plan is $%.2f.\n', TC)
TCp = sum(sum(Cp.*Fp));
TCi = sum(sum(Ci.*Fi));
TCf = sum(sum(Cf.*Ff));
vdisp('TCp,TCi,TCf,TC')