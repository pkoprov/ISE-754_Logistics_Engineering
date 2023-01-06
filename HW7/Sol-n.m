DCXY = uszip5('XY', mand([20548, 26149,36317],uszip5('Code5')));
CXY = uszip5('XY', mand([30669, 38339, 30732, 23830, 23154],uszip5('Code5')));
dem = [40 55 35 70 25];
D = dists(DCXY,CXY,'mi');
F = full(sparse(argmin(D),1:size(D,2), dem));
mdisp(F)
%%
sup = [60 90 80];
F = trans(D, sup,dem)
%% 
IJC = lev2list(D);
IJCU = [IJC repmat(30,size(IJC,1),1)];
s = [sup -dem];
lp = mcnf2lp(IJCU,s);
[x,TC,XFlg,out] = lplog(lp{:});
f = lp2mcnf(x,IJCU,s);
F = reshape(f,size(D));
mdisp(F)