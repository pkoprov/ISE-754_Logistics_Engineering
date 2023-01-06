%% Exam1
%% Problem 1
clear
disp('Problem 1');
%% Define initial conditions
%% Geolocate aggreagate demand points
NC = uszip5(strcmp('NC', uszip5('ST')));
OC = uscounty(mand('NC', uscounty('ST'), 'Orange', uscounty('Name')));
s = uscenblkgrp(uscenblkgrp('SCfips') == OC.SCfips);
OCxy = [-79.257497 36.243968; % Coordinates of Orange county limit
    -78.951653 36.239008;
    -79.018391 35.861254;
    -79.249062 35.878076];

makemap(OCxy)
% pplot(NC.XY, 'r.')
OC_ctr = pplot(OC.XY, 'bv')
CBG = pplot(s.XY, 'go')
BRD = pplot([OCxy; OCxy(1,:)], 'bs-');
%% Define the demand and distances
dem = s.Pop*4/2000;

dafh = @(XY1,a1,XY2,a2) max(1.2*dists(XY1,XY2,'mi'),...
   0.675*max(sqrt(a1(:)),sqrt(a2(:)')));

Da = dafh(s.XY,s.LandArea,s.XY,s.LandArea);  % Area-adjusted distances
%% Build constatnt and variable cost matrices
ppiTL = 136.3; % Jan 2020 (P)
r = 2*(ppiTL/102.7);
C = r*Da.*dem;
V = 8e4/2000; % maximum capacity for the facility
k = repmat(1000,1,size(C,1)); % fixed cost = 50
%% Create MILP model for CFL
clear mp
mp = Milp('CFL')
[n m] = size(C);
mp.addobj('min',k,C)
for i = 1:n
    mp.addcstr({V,{i}},'>=',{dem,{i,':'}})    
end
for j = 1:m
   mp.addcstr(0,{':',j},'=',1)
   
end
mp.addub(1,1)
mp.addctype('B','C')
mp.Model
%% Solve using Gurobi
clear model param
model = mp.milp2gb;
params.outputflag = 1;
result = gurobi(model, params);
x = result.x;
x = mp.namesolution(x);
TC = result.objval;
idxNF = find(round(x.k));
nNF = sum(x.k);
NF = pplot(s.XY(idxNF,:), 'kx')
legend([OC_ctr CBG BRD NF], ["Orange County center"...
    "Census Block Group" "County border" "New Facilities"])
fprintf('The optimal number of centers with 0 fixed cost is %d.\n', nNF);
fprintf('They should be located at:\n');
lonlat2city(s.XY(idxNF,:), uscity)

%% Problem 2
clear
disp('Problem 2');
%% Define initial conditions
f = [ 720 120 480];
wt = [16 150 40];
cu = [2 6 12];
uc = [500 200 100];
sv = [50 180 80];
Sup = uscity(mand('Chicago',uscity('Name'),'IL',uscity('ST')));
Cust = uscity(mand('Austin',uscity('Name'),'TX',uscity('ST')));
makemap([Sup.XY; Cust.XY]);
pplot(Sup.XY, 'ro'), pplot(Sup.XY, Sup.Name);
pplot(Cust.XY, 'bs'), pplot(Cust.XY, Cust.Name);
%% Define shipment and truck structures
ppiTL = 136.3; % Jan 2020 (P)
ppiLTL = 193.6; % Jan 2020 (P)
tr = struct('r',2*(ppiTL/102.7),'Kwt',25,'Kcu',2750);
hobs = (uc-sv)./uc; %obsolescence rate
h = 0.11+hobs;
sh = vec2struct('f', f.* wt/2000, 's', wt./cu, 'd',...
    1.2*dists(Sup.XY, Cust.XY, 'mi'), 'v', uc, 'h', h, 'a', 1);
sdisp(sh);
%% find the optimal TLC for separate shipping
[~,q,isLTL] = minTLC(sh,tr,ppiLTL);
[TLC,TC,IC] = totlogcost(q,transcharge(q,sh,tr),sh);
days = 365.25*q./[sh.f];
%% find the optimal TLC for aggregate shipping
ash = aggshmt(sh);
[~,qa,isLTLa] = minTLC(ash,tr,ppiLTL);
[TLCa,TCa,ICa] = totlogcost(qa,transcharge(qa,ash,tr),ash);
daysa = 365.25*qa./[ash.f];
%% display results
M = [ TLC sum(TLC); TC sum(TC); IC sum(IC); q NaN; days NaN; isLTL NaN]';
Ma = [TLCa TCa ICa qa daysa isLTLa'];
mdisp([M; Ma],{'1','2','3', 'Sum', 'Aggregate'},...
   {'TLC', 'TC', 'IC','q', 'days', 'isLTL'},'sh');

fprintf('\nThe optimal way to ship products is via aggregate TL.\n')
fprintf(['Aggregate shipment should be made every %.0f days'...
    ' with a size of %.2f tons.\n'], daysa, qa);
