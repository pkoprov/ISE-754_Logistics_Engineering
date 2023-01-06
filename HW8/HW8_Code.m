%% Q1b
clear
C = [
    0 8 6 9 1 5
    3 0 1 5 4 2
    9 2 0 3 1 1
    8 2 1 0 10 6
    6 7 10 1 0 10
    6 2 5 2 1 0]
loc = [1 5 2 6 4 3 1];
[loc, TC] = tsp2opt(loc,C);

fprintf(['The final location sequence determined after applying the'...
    'twoopt improvement procedure is %s, with total cost = %d.\n'],...
    num2str(loc), TC);
%% Q2b
clear
T = [
    0 2 2 2 1 2
    2 0 3 2 3 3
    2 3 0 3 2 2
    2 2 3 0 3 1
    1 3 2 3 0 3
    2 3 2 1 3 0];
mdisp(T)
sh = vec2struct('b',1,'e',[4 2 5 3 6]);
sh = vec2struct(sh,'tU',0, 'temin',[9 9 15 18 21],'temax',[18 12 18 21 24]);
tr = struct('b',1,'e',1,'tbmin',6,'tbmax',24,'temin',6,'temax',24);
sdisp(sh)
[TC,Xflg,out] = rteTC([1 2 3 4 5 1 2 3 4 5],sh,T,tr);
TC,Xflg
sdisp(out,false);
fprintf(['The minimum total time span needed to complete all '...
    'deliveries and return to the depot is %d hours.\n'], TC)
%% Q3b
clear
r123 = [2 3 2 1 3 1];
D = [
    0 180 320 100 100 40
    180 0 140 80 240 140
    320 140 0 220 300 280
    100 80 220 0 240 60
    180 240 300 240 0 220
    40 140 280 60 220 0];
ppi = 125;
tr = struct('r',2,'Kwt',25,'Kcu',2750);
sh = vec2struct('f',[200 300 100],'s',[20 5 10],...
   'b',[6 3 2],'e',[5 1 4],'v',[20000 5000 10000],'a',1,'h',.3);
sh = vec2struct(sh,'d',diag(D([sh.b],[sh.e])));
sdisp(sh)

[TLC,q,isLTL] = minTLC(sh,tr,ppi,D,r123)
%% Q4
%% Add depot: Use end location of shipment 1 as depot for all shipments
clear
DC = table2array(readtable('HW8data.xlsx','Sheet', 1));
DC = flip(DC);
Cust = readtable('HW8data.xlsx','Sheet', 2);
XY_c = [Cust.Lon Cust.Lat];
XY = [DC; XY_c];
D = dists(XY,XY,'mi')*1.2;
q = Cust.Pkg;
s = 1;
maxtime = 7;

temin = [];
temax = [];
for i = 1:size(XY_c,1)
    if strcmp(Cust.T_W(i), 'M')
        temin = [temin 8];
        temax = [temax 12];
    elseif strcmp(Cust.T_W(i), 'A')
        temin = [temin 12];
        temax = [temax 17];
    else
        temin = [temin 17];
        temax = [temax 21];
    end
end

sh = vec2struct('b',1,'e',2:size(XY,1), 'q', q, 's', s);
sh = vec2struct(sh,'tU', 2/60,'temin',temin,'temax',temax);

tr = struct('b',1,'e',1, 'Kcu',99999, 'Kwt', 35);

sdisp(sh)
%% Get road network
expansionAroundXY = 0.12;
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
%% Add connector roads from cities to road network
[IJD11,IJD12,IJD22] = addconnector(XY,XY2,IJD);
%% Convert road distances to travel times (needs to be after ADDCONNECTOR)
v.IR = 70;  % Rural Interstate highways average speed (mph)
v.IU = 50;  % Urban Interstate highways average speed (mph)
v.R = 45;   % Rural non-Interstate roads average speed (mph)
v.U = 20;   % Urban non-Interstate roads average speed (mph)
v.C = 15;   % Facility to road connector average speed (mph)

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
%% Construct & improve routes:
rTDh0 = @(rte) rteTC(rte,sh,T,tr);
rTDh = @(rte) myrteTC(rte,rTDh0,maxtime);
tic
IJS = pairwisesavings(rTDh,sh); toc
% sound(sin(1:3000));
tic
r = twoopt(savings(rTDh,sh,IJS),rTDh); toc
% sound(sin(1:3000));
%% Display route output structure
[TC,Xflg,out] = rTDh0(r);
Bars=[];
for i = 1:length(out), Bars=[Bars; [out(i).Depart(1) out(i).Depart(end)]]; end
gantt(num2cell(Bars,2))
fprintf('The number of vans needed for tomorrow’s deliveries is 11.\n')
%% Q5
%% Create Data
clear, close all
s = readtable('HW8data.xlsx','Sheet', 3);

Cust = uszip5('XY', 'Code5', mor([s.orig; s.dest], uszip5('Code5')));

b = [];
e = [];
for i=1:length(s.orig)
    b = [b; find(Cust.Code5 == s.orig(i))];
    e = [e; find(Cust.Code5 == s.dest(i))];
end

D = dists(Cust.XY, Cust.XY, 'mi')*1.2;

ppiLTL = 193.6; % Jan 2020 (P)
ppiTL = 136.3; % Jan 2020 (P)
r = 2*(ppiTL/102.7);

tr = struct('r',r,'Kwt',25,'Kcu',2750);
sh = vec2struct('idx',1:length(s.orig),'b',b, 'e',e,...
    'f',s.ud.*s.wt/2000,'s',s.wt./s.cu,'v',2000./s.wt.*s.uc,'a',0.5,'h',.3);
sh = vec2struct(sh,'d',diag(D([sh.b],[sh.e])));
[TLC1,q1,isLTL] = minTLC(sh,tr,ppiLTL);
sh = vec2struct(sh,'TLC1',TLC1,'q1',q1,'t1',q1./[sh.f],'isLTL',isLTL);
qmax = maxpayld(sh,tr);
sh = vec2struct(sh,'qmax',qmax);
sdisp(sh)
%% Independent shipments
plotshmt(sh,Cust.XY,[],tr,true)
%% Consolidated shipments
rTCh = @(rte) minTLC(sh,tr,ppiLTL,D,rte);
ph = @(rte) plotshmt(sh,Cust.XY,rte,tr);
IJS = pairwisesavings(rTCh,sh,minTLC(sh,tr));
%% Construct and improve routes
[rc,TLCc] = twoopt(savings(rTCh,sh,IJS,ph),rTCh,ph);
%% Make shipments not in routes into single-shipment routes
[rc,idx1,TLCc] = sh2rte(sh,rc,rTCh);
plotshmt(sh,Cust.XY,rc,tr)
%% Change in TLC from indep to consol:
100*(sum(TLCc) - sum(TLC1))/sum(TLC1)
%% Change in TLC for just multi-shipment routes
idxrte = find(cellfun(@length,rc) > 2);
idxsh = rte2idx(rc(idxrte));
idxsh = [idxsh{:}];
100*(sum(TLCc(idxrte)) - sum(TLC1(idxsh)))/sum(TLC1(idxsh))
%% Q6
clear, close all
load shmtNC30
tr = struct('r',2,'Kwt',25,'Kcu',2750);
idx = [1 3 26 5];
sh = vec2struct('idx',idx,'b',b(idx),'e',e(idx),'f',f(idx),...
    's',s(idx),'v',v(idx),'a',.5,'h',.3);
sh = vec2struct(sh,'d',diag(D([sh.b],[sh.e])));
sdisp(sh)
n = length(sh)
%%
rTCh = @(rte) minTLC(sh,tr,[],D,rte);
[TLC1,q1,isLTL] = minTLC(sh,tr,[])
%% Min incremental charge for all possible routes
R = perms(1:n)
R = sortrows(R,1:n)
C = zeros(size(R));
for i = 1:size(C,1)
   for j = 1:size(C,2)
      Rj = perms(R(i,1:j));  % Try all permutations to get optimal
      TC(j) = Inf;
      for k = 1:size(Rj,1)
         [~,TCj] = insertimprove(Rj(k,:),rTCh,sh);
         if TCj < TC(j), TC(j) = TCj; end
      end
   end
   C(i,:) = TC;
   TC = diff([0 TC]);
   C(i,:) = TC(invperm(R(i,:)));
end
mdisp(C,sum(R.*repmat(10.^[n-1:-1:0],size(R,1),1),2))
%% Equal charge allocation
TCc = min(sum(C,2))
c_equal = repmat(TCc/n,1,n)
pct_reduct = round(100*(1 - c_equal./TLC1))
%% Equal savings allocation
Sn = sum(TLC1) - TCc
c_eq_sav = TLC1 - Sn/n
pct_reduct = round(100*(1 - c_eq_sav./TLC1))
%% Exact Shapley allocation
c_Shap_exact = mean(C,1)
pct_reduct = round(100*(1 - c_Shap_exact./TLC1))
%% Pairwise approximate Shapely allocation
[~,S2] = pairwisesavings(rTCh,sh)
c_Shap_approx = TLC1 - (Sn/n + sum(S2)/(n-1) - sum(sum(S2))/(n*(n-1)))
pct_reduct = round(100*(1 - c_Shap_approx./TLC1))
%% Comparison
vdisp('TLC1,c_equal,c_eq_sav,c_Shap_exact,c_Shap_approx',true,true)