%%Q1
q= 25*70/2000;
uw=70;
uc=20;
s=uw/uc;
d=532;
ppiLTL=144.3;
rate = rateLTL(q,s,d,ppiLTL)
c_LTL=rate*q*d
%%
MC = 95.23; disc = 0; 
cl = [500 400 300 250 200 175 150 125 110 100 92.5 85 77.5 70 65 60 55 50];
cl_avg = [0.52 1.49 2.49 3.49 4.49 5.49 6.49 7.49 8.49 9.72 11.22 12.72 14.22 18.01 25.5 32.16 39.68 56.18];
class = cl(argmin(dists(s,cl_avg')));
qb = [0 0.25 0.5 1 2.5 5 10 15 20 Inf];
i = find(qb(1:end-1) <= q & q < qb(2:end));
ODi = 124.23;
ci = ODi*20*q;
ODiplus1 = 101.83;
qbi = qb(i+1);
ciplus1=ODiplus1*20*qbi;
c_tar = (1-disc)*max(MC, min(ci, ciplus1))
%% diff
c_tar/c_LTL
%% Q2
clear all
sh.f = 75;
sh.d= 625;
sh.a=1;
ppiTL = 123.4;
ppiLTL = 141.4;
sh.v = 11200;
sh.s = 12;
sh.h = 0.4;
tr.Kwt = 25; tr.Kcu = 2750; tr.ppiTL = ppiTL; tr.r = 2 * (ppiTL/102.7);
[TLC,q,isLTL] = minTLC(sh,tr,ppiLTL);
n = 7/365.25;
q1w = sh.f*n;
[c,isLTL,cTL,cLTL] = transcharge(q1w,sh,tr,ppiLTL) 
[TLC1w,TC1w,IC1w] = totlogcost(q1w,c,sh)
incr_in_TLC = TLC1w - TLC
%% Q3
clear all
uwt = [24 40];
ucu = [10 4];
udem = [1000 1200];
uval = [35 40];
a = 0.5;
th = [2 3];
sh = vec2struct('d', 500, 'f', udem.*uwt/2000, 's', uwt./ucu, 'a',a, 'v', uval./(uwt/2000), 'h',0.04+0.06+0.7./th);
tr.Kwt = 25; tr.Kcu = 2750;  tr.r = 2;
sdisp(sh);
[TLC,q] = minTLC(sh,tr);
vdisp('TLC,q')
%% aggship
assh = aggshmt(sh);
sdisp(assh);
[TLCa,qa] = minTLC(assh,tr);
vdisp('TLCa,qa')
%% diff
(TLCa-sum(TLC))/sum(TLC)
%% Q4
clear all
city = {'Asheville', 'Statesville', 'Winston-Salem','Greensboro', 'Durham', 'Raleigh','Wilmington'};
P = [50 150 190 220 270 295 420]';
shS = vec2struct('f', [400 120], 's', [100 50]./[25 5]);
shC = aggshmt(shS);
shC = vec2struct(shC,'f', shC.f*[20 25 30 15 10]/100);
sh = [shS shC];
sdisp(sh);
tr.Kwt = 25; tr.Kcu = 2750;
qmax = maxpayld(sh,tr);
n = [sh.f]./qmax;
vdisp('qmax,n');
x = minisumloc(P([2 7 1 3 4 5 6]),n,1)
city(argmin(dists(P,x)))
%% Q5
clear
city = {'Richmond', 'Canton', 'Malden','Tyler'};
st = {'CA', 'OH', 'MA', 'TX'};
zip = [72118, 55472, 15010, 88102, 87301, 73099];
ud = [1200, 3200, 2200, 1100, 1600, 1500];
cu = [2.7, 1.3, 2.7, 3.6];
wt = [8, 14, 4, 29];
BOM = [4, 3, 1, 4];
cuFG = 38.2; wtFG = 175;
%%
city2lonlat = @(city,st) uscity('XY', mand(city, uscity('Name'),st,uscity('ST')))
for i = 1:length(city)
    XYP(i,:) = city2lonlat(city{i}, st{i});
end
XYC = uszip5('XY', mand(zip, uszip5('Code5')));
%%
tr.Kwt = 25; tr.Kcu = 2750;
fout = ud*wtFG/2000;
qmaxFG = maxpayld(wtFG/cuFG,tr);
wout = fout/qmaxFG;
fin = (BOM*sum(ud)).*wt/2000;
qmax = maxpayld(wt./cu, tr);
win = fin./qmax;
xy = minisumloc([XYP;XYC], [win wout], 'mi');
cityst = uscity50k;
idx = lonlat2city(xy, cityst);
cityst.Name{idx}, cityst.ST{idx}
%% Q6
clear
fn = 'HW6data.xlsx'
inS = table2struct(readtable(fn, 'Sheet', 'Supplier'));
inC = table2struct(readtable(fn, 'Sheet', 'Customer'));
UD = xlsread(fn, 'Demand');
%%
city2lonlat = @(city,st) uscity('XY', mand(city, uscity('Name'),st,uscity('ST')));
DCstr = {'Durham', 'NC'};
DC = city2lonlat(DCstr{:});
SXY = uszip5('XY', mand([inS.zip], uszip5('Code5')));
CXY = uszip5('XY', mand([inC.zip], uszip5('Code5')));
XY = [SXY; DC; CXY];
makemap(XY);
hS = pplot(SXY, 'ro');
hDC = pplot(DC, 'gv');
hC = pplot(CXY, 'k.');
pplot(DC, DCstr(1));
legend([hS hDC hC], {'Suppliers', 'DC', 'Customers'});
%%
ppiTL = 135.1;
tr = struct('r', 2 * (ppiTL/102.7), 'Kwt', 25, 'Kcu', 2750);
F = [inS.wt]'.*UD/2000;
s = [inS.wt].*[inS.cu];
v = 2000*[inS.uc]./[inS.wt];
h_ob = ([inS.uc]-[inS.sv])./[inS.uc];
h = 0.05 +0.06+h_ob;
in.a0 = 0; in.aD = 0.5;
shS = vec2struct('f', sum(F,2), 's',s,'v',v, 'h',h,'d', dists(SXY, DC,'mi')*1.2);
shC = vec2struct(aggshmt(shS), 'f', sum(F,1), 'a', 0+in.aD, 'd', dists(DC,CXY,'mi')*1.2);
in.ppiLTL = 188.3;
%% no ccord
sh = [vec2struct(shS, 'a', in.a0+0.5) shC];
[NC.TLC,NC.q,NC.isLTL] = minTLC(sh,tr,in.ppiLTL);
NC.t = 365.25*[NC.q]./[sh.f];
vdisp('NC.TLC,NC.q,NC.t,NC.isLTL', true, true)
%% all cross-docked
sh = [vec2struct(shS, 'a', in.a0) shC];
qmax = maxpayld(sh,tr);
TLC0h = @(t) totlogcost([sh.f]*t,transcharge([sh.f]*t,sh,tr,in.ppiLTL),sh);
TLCh = @(t) iff(t>=1/365.25 & [sh.f]*t<=qmax, TLC0h(t), Inf);
tx0 = min(qmax./[sh.f]);
X.t = fminsearch(@(t) sum(TLCh(t)),tx0);
X.q = [sh.f]*X.t;
[~,X.isLTL] = transcharge(X.q,sh,tr,in.ppiLTL) ;
X.TLC = TLCh(X.t);
vdisp('X.TLC, X.q, X.isLTL', true, true)