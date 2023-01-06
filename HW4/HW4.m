%% HW 4
fprintf("HW4\nPavel Koprov\n")
%% Q1
disp("Question 1")
%%  17. Full vs. Sparse Matrices
disp("17. Full vs. Sparse Matrices")
A = [0 0 0 30; 10 0 20 0]
sparse(A)
full(A)
%%  18. Inserting and Extracting Array Elements
disp("18. Inserting and Extracting Array Elements")
A = sparse([1 2 2], [4 1 3], [30 10 20], 2, 4)
disp("Inserting Into an Existing Matrix")
B = [1:4; 5:8]
B(A ~= 0) = A(A ~= 0)
disp("Extracting From a Matrix")
B([1 2 2],[4 1 3])
b = diag(B([1 2 2],[4 1 3]))
%% Q4
disp("Question 4")
x = [10020 17112 27707 32606 48234 56123];
c = uszip5(mor(x ,uszip5('Code5')));
s = uscity(mor({'Gainesville', 'Warren'}, uscity('Name'))...
    & mor({'FL', 'OH'}, uscity('ST')));
P = [s.XY(1:2,:); c.XY];
makemap(P)
pplot(c.XY, 'r.')
pplot(s.XY(1:2,:),'g*'), pplot(s.XY(1:2,:),s.Name(1:2,:))
fo = [10 40 35 15 30 25];
fi = [1600 250] *sum(fo)/2000 ;
ro = 0.06;
ri = [0.012 0.025];
wo = ro * fo;
wi = fi .* ri;
w = [wi wo];
%% a)
disp("a)")
rng(1234)
[X, TC] = minisumloc(P,w,'mi');
pplot(X, 'bv')
cty = uscity50k(argmin(dists(X,uscity50k('XY'),'mi')));
pplot(cty.XY, 'kx'), pplot(cty.XY, cty.Name);
fprintf(['The nearest 50K city to optimal location of the widget factory'...
    'is %s, %s.\n'], char(cty.Name), char(cty.ST));
%% b)
disp("b)")
f12 = sum(fo(1:4));
f13 = sum(fo(5:6));
w12 = f12*0.04;
w13 = f13*0.04;
V = full(sparse([1 1], [2 3], [w12 w13], 3,3));
W = full(sparse([1 1 2 2 2 2 3 3], [1:8], [wi(1:2) wo(1:6)],3,8));

rng(1234)
[X1,TC1] = minisumloc(P,W,'mi',V)

for i = 1:length(X1)
    cty = uscity50k(argmin(dists(X1(i,:),uscity50k('XY'),'mi')));
    fprintf(['The nearest 50K city to optimal location of the'...
        'facility %d is %s, %s.\n'], i, char(cty.Name), char(cty.ST));
end
pplot(cty.XY,'gv')
pplot(cty.XY, cty.Name)
%% c)
disp("c)")
dTC = (TC-TC1)/TC;
fprintf(['The percentage change in total cost associated with using'...
    ' DCs is %.3f%%.\n'], (dTC*100))

%% Q5
disp("Question 5.b)")
P = [50 150 190 220 270 295 420]';
m = size(P,1);
w = ones(1,m);
X = [60 125 130]';
n = size(X,1);
p = 1;
rng(1234)
[xloc,TC] = ala(X,w,P,p);
disp('The optimal locations of 3 facilities are')
xloc
%% Q6
disp("Question 6)")
DC = uscity(strcmp('Roanoke Rapids', uscity('Name')));
s = uscity;
ind = find(strcmp('Charlottesville',s.Name));
CharXY = s.XY(ind,:);
CharName = s.Name(ind,:);
s = uscity(mor({'NC','SC','VA'},uscity('ST')) & uscity('Pop')>5000 &...
    s.XY(:,2) < CharXY(2));

makemap(s.XY)
pplot(s.XY, 'r.')
pplot(DC.XY, 'r*')
pplot(CharXY, 'g*')
pplot(CharXY, CharName)
TC_Cary = 5.7e+6;
fo = s.Pop;
P = s.XY;
r = TC_Cary/sum(fo'.*dists(DC.XY,P,'mi'));
w = fo*r;
%% (a) What would be the maximum expected reduction in annual
% outbound transportation costs if the DC could be re-located to 
% any other location?
rng(1234)
[Xa, TCa] = minisumloc(s.XY,w','mi');
pplot(Xa, 'bv');
pplot(Xa, {'a)'});
dTC = TC_Cary - TCa;
fprintf(['(a) The maximum expected reduction in annual outbound'...
    '  transportation costs is %.2f$.\n'], dTC);
%% (b) What would be the maximum expected reduction in annual
% outbound transportation costs if three DCs could be located
% anywhere and the existing DC would be closed, where each
% store would be supplied by one of the DCs?
n = 3;
rng(1234)
[Xb,TCb] = ala(randX(P,n),w,P,'mi');
% pplot(Xb, 'bv');
% pplot(Xb, {'b)'});
dTC = TC_Cary - TCb;
fprintf(['(b) The maximum expected reduction in annual outbound'...
    '  transportation costs is %.2f$.\n'], dTC);
%% (c) What would be the maximum expected reduction in annual outbound
% transportation costs if two DCs could be located anywhere and the
% existing DC would be closed, where each store would be supplied 
% by one of the DCs? Also compare with the cost obtained in (b).
n = 2;
rng(1234)
[Xc,TCc] = ala(randX(P,n),w,P,'mi');
dTC = TC_Cary - TCc;
fprintf(['(c) The maximum expected reduction in annual outbound'...
    '  transportation costs is %.2f$.\n'], dTC);
dTCb = TCb-TCc;
fprintf([' Compared to (b), (c) is %.2f$ more expensive.\n'], abs(dTCb));
%% (d) What would be the change from (c) if only one new DC could be
% located anywhere and the existing DC in Roanoke Rapids remained open?
rng(1234)
loc_h = @(W, X0) [DC.XY; minisumloc(P,W(2,:),'mi',[],X0(2,:))];
[Xd,TCd] = ala(randX(P,n),w,P,'mi',loc_h);
dTC = TCc - TCd;
fprintf(['(d) The change from (c) is %.2f$.\n'], dTC)