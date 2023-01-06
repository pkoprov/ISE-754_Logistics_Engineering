%% Q1
disp('Question 1')
P = [85 40; 75 40; 75 85; 45 35; 40 55; 35 55; 30 85; 55 90];
w = [11 46 36 20 32 23  20 18];
d1h = @(x, P) sum(abs(x - P), 2);
TDh = @(x) sum(w(:) .* d1h(x, P));
xy = fminsearch(TDh, mean(P))
%% Q2
disp('Question 2')
fout = [10 38 20 46];
rout = 0.8;
wout = rout * fout;
BOM = [2 0.55 1.2 0.35];
fin = BOM * sum(fout);
rin = [0.08 0.05 0.15 0.03];
win = fin .* rin;
w = [win wout];
%% a)
disp('a)')
P = [270 150 420 50 50 190 220 295]';
d2h = @(x, P) sum(abs(x - P), 2);
TCh = @(x) sum(w(:) .* d2h(x, P));
x = fminsearch(TCh, mean(P));
fprintf('Optimal location is on the mile marker %d.\n', round(x))
%% b)
disp('b)')
disp(['The product is '...
    iff(sum(win) < sum(wout), 'Weight-gaining', 'Weight-loosing')])
%% Q3
disp('Question 3')
R = [2 1; 10 7];
rng(2073);
P = R(1, :) + diff(R) .* rand(20, size(R, 2));
d3h = @(x, P) sqrt(sum((x - P) .^2, 2));
plot(P(:,1), P(:,2), '.r', 'Displayname', 'EFs'), shg
if ~exist('iff')
    fid = fopen('iff.m', 'wt');
    fprintf(fid, 'function x=iff(a,b,c)\n if a, x=b;else x=c;end;\n');
    fclose(fid);
end
%% (a) Determine the location that minimizes the sum of the squared-Euclidean distances.
disp('a)')
TCh = @(x) sum(d3h(x, P) .^2);
x = mean(P,1);
[xy, TC] = fminsearch(TCh, x)
hold on
plot(xy(1), xy(2), 'g*', 'Displayname', 'Q3.a'), shg
%% (b) Determine the location that minimizes the sum of the Euclidean distances.
disp('b)')
TCh = @(x) sum(d3h(x, P));
x = mean(P,1);
[xy, TC] = fminsearch(TCh, x)
hold on
plot(xy(1), xy(2), 'bd', 'Displayname', 'Q3.b'), shg
%% (c) Determine the location that minimizes the maximum Euclidean distance.
disp('c)')
TCh = @(x) max(d3h(x, P));
x = mean(P,1);
[xy, TC] = fminsearch(TCh, mean (P,1))
hold on
plot(xy(1), xy(2), 'c*', 'Displayname', 'Q3.c'), shg
%% (d) Determine the location that minimizes the sum of the Euclidean distances, with the
% restriction that the location cannot be within the rectangle defined by (5,3) and (7,5). The
% straight line from the location to any point can cross the rectangle (e.g., if the location
% restriction is due to zoning, then travel across the zone can still occur).
disp('d)')
Zone = [5 3; 7 5];
rectangle('Position', [Zone(1,:) diff(Zone, [], 1)])
inZone = @(x, Z) x(1) >= Z(1,1) && x(1) <= Z(2,1) && x(2) >= Z(1,2) && x(2) <= Z(2,2);
TCh = @(x) iff(inZone(x, Zone), Inf, sum(d3h(x, P)));
n = [2 1; 2 7; 10 7; 10 1];
TC = Inf;
for i = 1:4
    [xi, TCi] = fminsearch(TCh, n(i,:));
    if TCi < TC, x = xi; TC = TCi;
    end
end
[x TC]
hold on
plot(x(1), x(2), 'bo', 'Displayname', 'Q3.d'), shg
%% (e) Determine the location that maximizes the minimum Euclidean distance, with the
% restriction that the location must be within the rectangle defined by (2,1) and (10,7).
disp('e)')
TCh = @(x) iff(inZone(x, R), -min(d3h(x, P)), Inf);
n = [2 1; 2 7; 10 7; 10 1];
TC = Inf;
for i = 1:4
    [xi, TCi] = fminsearch(TCh, n(i,:));
    if TCi < TC, x = xi; TC = TCi;
    end
end
[x TC]
hold on
plot(x(1), x(2), 'rd', 'Displayname', 'Q3.e'), shg
legend('show')
%% Q4
close all
a = [0 30];
w = [1 2];
d4h = @(x) abs(a - x);
k = 3
TCh = @(x) sum(w.*d4h(x).^k);
[x TC] = fminsearch(TCh, 0)
hold on
for i = 0:30
    plot(i, TCh(i), 'r*'), shg
end
plot (x, TCh(x), 'gd', 'Displayname', 'Min TC')
%%
close all
hold on
for i = 1:200
    k = i;
    TCh = @(x) sum(w.*d4h(x).^k);
    [x TC] = fminsearch(TCh, 0);
    plot(k,x, 'r*')
end
disp(['X is equal to ' num2str(x) ' while k is ' num2str(k)])