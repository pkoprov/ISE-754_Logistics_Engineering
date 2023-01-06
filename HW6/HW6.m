%% HW 6
%% Question 2
disp("Question 2")
clear
ppiTL = 133.2;
% ppiLTL = 141.4;
sh.d = 250;
sh.s = 4;
tr.Kwt = 25; tr.Kcu = 2750; tr.ppiTL = ppiTL;
sh.h = 41.74+0.11;
sh.f = 100;
sh.v = 1000;
sh.a = 1;
sdisp(sh)
%% Shipping at a constant rate without contraints on shipping intervals
tr.r = 2 * (ppiTL/102.7);
qFTL = maxpayld(sh,tr);
nFTL = sh.f/qFTL;
days = 365.25*qFTL./sh.f;
[c,isLTL,cTL] =  transcharge(qFTL,sh,tr);
[TLC_FTL,TC_FTL,IC_FTL] = totlogcost(qFTL,c,sh);
disp('TLC if transporting at a constant TL/year rate')
vdisp('TLC_FTL,qFTL,TC_FTL,IC_FTL,isLTL,nFTL,days')
%% If sending in LTL every week
tmax = 7/365.25/2;
nmin = 1/tmax;
qw = sh.f/nmin;
[c,isLTL,cTL] =  transcharge(qw,sh,tr);
[TLCw,TCw,ICw] = totlogcost(qw,c,sh);
disp('TLC if transporting once a week in LTL')
vdisp('TLCw,qw,TCw,ICw,isLTL, nmin')
%% change in TLC
dTLC = TLC_FTL - TLCw;
fprintf('The change in TLC is $%.2f.\n', dTLC)
%% Question 3
fprintf("\nQuestion 3")
clear
sh = vec2struct('d', 500,'f', 1000*24/2000, 's', 24/10,...
    'a', 1, 'h', (70/2+4+6)/100);
sh(2) = vec2struct(sh(1), 's', 40/4, 'h', (70/3+4+6)/100, 'f', ...
    1200*40/2000);
sh(1).v = 35*1000/sh(1).f; sh(2).v = 40*1200/sh(2).f;
sdisp(sh)

tr.Kwt = 25; tr.Kcu = 2750; tr.r = 2;
%% Separate shipping
qmax = maxpayld(sh,tr);
[TLC,TC,IC] = totlogcost(qmax,transcharge(qmax,sh,tr),sh);
t = 365.25*qmax./[sh.f];
fprintf('\nTLC if transporting A and B separately')
vdisp('TLC,qmax,TC,IC, t', true, true);

%% Aggregate Shipment
ash = aggshmt(sh);
fprintf('\nAggregate Shipment');
sdisp([sh ash],[],'sh12,agg3');
[TLCa,qa] = minTLC(ash,tr);

M = [ TLC sum(TLC); maxpayld(sh,tr) NaN; t NaN]';
Ma = [TLCa maxpayld(ash,tr) 365.25*qa./ash.f'];
fprintf('\nTLC if transporting A and B together')
mdisp([M; Ma],{'A','B','A+B','Aggregate'},...
   {'TLC','qmax', 't'},'sh')
%% change in TLC
dTLC = sum(TLC) - TLCa;
fprintf('The change in TLC is $%.2f.\n', dTLC)
%% Question 5
fprintf("\nQuestion 5")
clear
s = 175/38.2;
fC = [1200 3200 2200 1100 1600 1500];
shC = vec2struct('f',fC*175/2000,'s',s);
fS = [4 3 1 4]*sum(fC);
uwt = [8 14 4 29];
ucu = [2.7 1.3 2.7 3.6];
shS = vec2struct('f',fS.*uwt/2000,'s',uwt./ucu);
sh = [shS shC];
sdisp (sh)

city2lonlat = @(city,st) ...
   uscity('XY',mand(city,uscity('Name'),st,uscity('ST')));
ZIP2lonlat = @(ZIP) uszip5('XY',mand(ZIP,uszip5('Code5')));
Sxy = city2lonlat(...
    ["Richmond","Canton","Malden", "Tyler"]', ["CA", "OH", "MA", "TX"]');
Cxy = ZIP2lonlat([72118, 55472, 15010, 88102, 87301, 73099]);
%% FTL: Locate DC to min transport cost assuming FTL used for all transport
tr = struct('r',2,'Kwt',25,'Kcu',2750);
qmax = maxpayld(sh,tr);
n = ([sh.f]./qmax);
w = n*tr.r;
DC = minisumloc([Sxy; Cxy],w,'mi');
[idx,dist,drt,dstr] = lonlat2city(DC);
fprintf('Plant should be located in %s.\n', dstr(end-8:end))
%% Plot the facilities
makemap([Sxy; Cxy])
pplot(Sxy, 'rv'), pplot(Sxy, {'Richmond','Canton','Malden', 'Tyler'});
pplot(Cxy, 'ro'), pplot(Cxy,{'72118', '55472', '15010', '88102', '87301', '73099'})
pplot(DC, 'bx', 'MarkerSize',12), pplot(DC, uscity50k('Name', idx))
%% Question 6
fprintf("\nQuestion 6")
clear
DC = uscity(mand('Durham', uscity('Name'), 'NC', uscity('ST')));
Sup = readtable('HW6data.xlsx','Sheet', 1);
Dem = readtable('HW6data.xlsx','Sheet', 2);
Cust = readtable('HW6data.xlsx','Sheet', 3);
Cust = uszip5(mand(table2array(Cust),uszip5('Code5')));
SupZip = uszip5(mand(Sup.zip,uszip5('Code5')));
Dem = table2array(Dem);

uwt = Sup.wt;                   %unit weight of products
ucu = Sup.cu;                   %unit weight of products
hobs = (Sup.uc-Sup.sv)./Sup.uc; %obsolescence rate
h = 0.11+hobs;                  %inventory carrying rate
Sdst = dists(DC.XY, SupZip.XY, 'mi')'*1.2; % distances from Sup to DC
fS = sum(Dem,2);                %demand for every product
vS = (Sup.uc.*fS)./(fS.*uwt/2000);  %cost of product per ton
shS = vec2struct('f', fS.*uwt/2000, 'd', Sdst, 's', uwt./ucu,...
    'v', vS, 'h', h, 'a', 0);

fC = sum(Dem.*uwt/2000,1)';
% vdisp('sum(fS.*uwt/2000), sum(fC)')      %checking the correctness
Cdst = dists(DC.XY, Cust.XY, 'mi')'*1.2; %distances from DC to Cust
sagg = sum(Dem.*uwt,1)./sum(Dem.*ucu,1); %aggregate density of outbound shipment
vC = sum(Sup.uc.*Dem,1)'./fC;
shC = vec2struct('f', fC, 'd', Cdst, 's', sagg',...
    'v', vC, 'h', 0.3, 'a', 1/2);

sh = [shS shC];
% sdisp(sh)

ppiTL = 131.4; % Jan 2018 (P)
ppiLTL = 179.4; % Jan 2018 (P)
tr = struct('r',2*(ppiTL/102.7),'Kwt',25,'Kcu',2750);
%% Plot
makemap([DC.XY; Cust.XY; SupZip.XY])
hDC = pplot(DC.XY, 'bo');
hC = pplot(Cust.XY, 'r.');
hS = pplot(SupZip.XY, 'g.');
legend([hS hDC hC],{'Suppliers','DC','Customers'});
%% 2. Use a DC with No Coordination
%% Determine optimal shipment sizes
qmax = maxpayld(sh,tr);
[TLC2,q2,isLTL2] = minTLC(sh,tr,ppiLTL);
n2 = [sh.f]./q2;
t2 = 365.25./n2;
vdisp('TLC2,q2,qmax,n2,t2,isLTL2',true,true)
%% 3: Use a DC with Perfect Cross-Docking
%% Determine optimal common shipment interval
TLC30h = @(t) totlogcost([sh.f]*t,transcharge([sh.f]*t,sh,tr,ppiLTL),sh);
TLC3h = @(t) iff([sh.f]*t <= qmax, TLC30h(t), Inf);
tx0 = min(qmax./[sh.f]);
t3 = fminsearch(@(t) sum(TLC3h(t)),tx0)
%% Report results
q3 = [sh.f]*t3;
[~,isLTL3] = transcharge(q3,sh,tr,ppiLTL);
TLC3 = TLC3h(t3);
n3 = [sh.f]./q3;
t3 = 365.25./n3;
vdisp('TLC3,q3,qmax,n3,t3,isLTL3',true,true)
%% change in TLC
dTLC = sum(TLC3) - sum(TLC2);
fprintf('\nThe change in TLC is $%.2f.\n', dTLC)
disp('All products should be stocked with no coordination.')