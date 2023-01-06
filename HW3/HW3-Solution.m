%% 2a)
s = uscity10k(mor({'NC', 'SC'}, uscity10k('ST')));
TCh = @(xy) sum(dists(xy,s.XY,'mi') * s.Pop);
xy = fminsearch(TCh,mean(s.XY))
xy = minisumloc(s.XY, s.Pop', 'mi')
%% b)
lonlat2city(xy, uscity10k)
d = dists(xy,s.XY,'mi');
s.Name(argmin(d))
%% c)
xy1 = minisumloc(s.XY, (-s.Pop)', 'mi')
lonlat2city(xy, s)
%
s.Name(argmax(d))
%% d
idx = argsort(d)
s.Name(idx(1:4))
%% e)
makemap(s.XY)
pplot(s.XY(idx(1:4),:), s.Name(idx(1:4)))
pplot(s.XY(idx(1:4),:), 'k*')
%% f)
dists(xy,s.XY(argmax(s.Pop),:), 'mi')
lonlat2city(s.XY(argmax(s.Pop),:))
%% g)
pool = s.Name(s.Pop>=50000)
pool(argmin(dists(xy, s.XY(s.Pop>=50000,:), 'mi')))
%% h)
TotalPop = sum(s.Pop)
SouthPop = sum(s.Pop(s.XY(:,2)<xy(2)))
perc = SouthPop/TotalPop*100
%% j)
range100Pop = sum(s.Pop(d<=100))
%% k)
s50 = uscity10k(uscity10k('Pop')>=50000);
d50 = dists(xy,s50.XY,'mi');
idx = argsort(d50);
makemap(s50.XY(idx(1:5),:))
pplot(s50.XY(idx(1:5),:), s50.Name(idx(1:5)))
pplot(s50.XY(idx(1:5),:), 'k*')
%% l)
s(Name)
idx = mand({'Ral', 'charl', 'north char'}, s.Name,{'NC', 'NC', 'SC'}, s.ST);
D = dists(s.XY(idx,:),s.XY, 'mi');
sum(sparse(argmin(D,1),1:length(s.Pop),s.Pop),2)
%% 3)
W = [20, 30, 24];
[Name, XY, Pop] = uscity('Name', 'XY', 'Pop', mand({'Detroit', 'Gainesville', 'Memphis'}, uscity('Name'), {'MI','FL','TN'}, uscity('ST')))
%% a)
xy = minisumloc(XY, W, 'mi')
%% b)
lonlat2city(xy)
%% c)
lonlat2city(xy, uscity10k)
%% d)
d = dists(XY,XY,'mi')
g = mean([1057/d(1,2), 754/d(1,3), 719/d(2,3)])
%% e)
TC = sum(W.*dists(xy, XY,'mi')*3*g)
%% f)
carXY = uscity('XY', mand('Cary', uscity('Name'), 'NC', uscity('ST')))
TC_Cary = sum(W.*dists(carXY, XY,'mi')*3*g)
diff = TC_Cary-TC
vdisp('TC, TC_Cary, diff')