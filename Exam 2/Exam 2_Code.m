%% Exam_2
%% Pavel 'Pasha' Koprov
%% Q1
clear
XY = uscity('XY', mand({'Amarillo' 'Anderson'}, uscity('Name'),...
    {'TX' 'SC'}, uscity('ST')));

cu = [6 3];         % unit cubic volume ft3
wt = [12 96];       % unit weight lb

T = 13;
D = [64 64 56 31 166 62 53 96 126 143 39 45 55;
    125 101 40 26 89 104 87 22 57 81 297 163 113]'; % unit demand 
d = dists(XY(1,:),XY(2,:),'mi')*1.2;
ppiTL = 136.3; % Jan 2020 (P)
r = 2*(ppiTL/102.7);
tr = struct('r',r,'Kwt',25,'Kcu',2750);
s = wt./cu
Ct = tr.r*d;         % transportation cost per truck per ton
qmax = maxpayld(s,tr)
Q = qmax
Cp = [120 80]      % unit cost to fabricate

yinit = [77; 100];   % initial storage 
yfinal = yinit;     % final storage 
hobs = [0.2 0.1]/13;
h = hobs + 0.11;
Ci = cumsum(Cp+Ct./Q).*h  % inventory cost for stage each product
T = 13;             % number of periods of production = 13
G = 2;              % number of products produced = 2
%% Create MILP model
Cp = repmat(Cp',1,T);    % create T x G array
Ci = repmat(Ci',1,T+1);
Ci(:,1) = 0 ; % Intital inventory cost already accounted for last period 
Ct = repmat(Ct',1,T);
mp = Milp('TPlan');
mp.addobj('min',Cp,Ci,Ct);  % Objective
for t = 1:T  % Flow balance constraints
    for m = 1:G-1
      mp.addcstr({[1 -1],{[m m+1],t}},{[1 -1],{m,[t t+1]}},'=',0)
   end
   mp.addcstr({G,t},{[1 -1],{G,[t t+1]}},'=',D(t));
   for g = 1:G
   mp.addcstr({g,t},0,'<=',{Q(g),{g,t}});
   end
end
mp.addlb(0,[yinit zeros(G,T-1) yfinal], 0)  % Lower bounds
mp.addub(Inf,[yinit inf(G,T-1) yfinal], Inf)    % Upper bounds
%% Solve using Gurobi
lp = mp.milp2lp;
[x,TC,exitflag,output] = linprog(lp{:});
TC, output
x = mp.namesolution(x)
%% Report results
Fp = x.Cp; mdisp(Fp)
Fi = x.Ci; mdisp(Fi)
Ft = x.Ct; mdisp(Ft)
mdisp(D')

%% Q2
clear, close all
df = table2struct(readtable('Exam2DataF20.xlsx'));

XY = [[df.Longitude]' [df.Latitude]'];
q = [df(2:end).Weight]'/2000;
s = [df(2:end).Density]';
tL = 20/60;
tU = 5/60;

sh = vec2struct('b',1,'e',[df(2:end).Customer]', 'q', q, 's', s);
tr = struct('b',1,'e',1,'tbmin',7,'temax',17,'Kwt',25,'Kcu',2750,...
    'maxTC', 10);

i = find([sh.q]'*2000./[sh.s]'/tr.Kcu > 1) % which shipment is above truck cubic capacity
srpls = sh(i).q*2000/sh(i).s/tr.Kcu - 1 % fraction of the surpluss

sh(end+1) = sh(i); % add additional shipment
sh(end).q = srpls*sh(i).q;

sh(i).q = sh(i).q - sh(end).q; % subtract additional shipment from overcubic shipment

sh = vec2struct(sh,'tU',[sh.q]'*tU, 'temin',7,'temax',17);
sdisp(sh)
%% Get road network
expansionAroundXY = 0.1;
[XY2,IJD,isXY,isIJD] = subgraph(usrdnode('XY'),...
   isinrect(usrdnode('XY'),boundrect(XY,expansionAroundXY)),...
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
%% Add connector roads from customers to road network
[IJD11,IJD12,IJD22] = addconnector(XY,XY2,IJD);
h = [h pplot(IJD12,[XY; XY2],'b-','DisplayName','Connector Roads')];
h = [h pplot(XY(2:end,:),'r.','DisplayName','Customers')];
h = [h pplot(XY(1,:),'g.','DisplayName','DC')];
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
%% Shortest time routes
n = size(XY,1);
[T,P] = dijk(list2adj([IJT12; IJT22]),1:n,1:n);
T = T+5/60;
%% Construct & improve routes:
rTCh = @(rte) rteTC(rte,sh,T,tr);
tic
IJS = pairwisesavings(rTCh,sh); toc
tic
r = twoopt(savings(rTCh,sh,IJS),rTCh); toc
%% add any single-shipment routes
[r,~,Time] = sh2rte(sh,r,rTCh);
%% Plot routes
plotshmt(sh,XY,r,tr)
%% Display route output structure
[TC,Xflg,out] = rTCh(r);
for i = 1:length(out), sdisp(out(i),false,i), end
%% Display Gantt chort of route spans
b= arrayfun(@(x) (x.Start(1)),out); b = b(:);
e= arrayfun(@(x) (x.Depart(end)),out); e = e(:);
figure
gantt([b e])
%% Route time and delivery cubic ft
for i = 1:length(r)
    idx = r{i}(isorigin(r{i}));
    Maxload(i) = sum([sh(idx).q]'*2000./[sh(idx).s]');
end
vdisp('Time, Maxload')
%% number of trucks
m = length(Time);
M = 1:m;
V = 10;
v = Time;
mp = Milp('# of Trucks');
mp.addobj('min',ones(1,m),zeros(m));
for i = M
   mp.addcstr({V,{i}},'>=',{v+tL,{i,':'}})   % V*yi >= sum_j(vj*xij)
end
for j = M
   mp.addcstr(0,{':',j},'=',1)            % sum_i(xij) = 1
end
mp.addctype('B','B')
%% Use INTLINPROG to solve
ilp = mp.milp2ilp;
x = intlinprog(ilp{:});
x = mp.namesolution(x);
B = arrayfun(@(i) find(x.arg2(i,:)),find(x.arg1),'UniformOutput',false);
B{:}
fprintf('Number of required trucks is %d.\n', length(B))
%% Chart for Trucks;

b=[]
for i=1:length(B)
    b = [b; 7 7+Time(B{i}(1))]
    for j=2:length(B{i})
        b = [b; b(end)+tL b(end)+tL+Time(B{i}(j))]
    end
end

figure
gantt([b])