%% HW3
%% 1. Duplicate the results of Sections 12–15 of Basic Concepts in Matlab.
disp('Q1. Duplicate the results of Sections 12–15 of Basic Concepts in Matlab.')
%% 12 Cell Arrays, Structures, and N-D Arrays
disp(" 12 Cell Arrays, Structures, and N-D Arrays")
%% Cell Arrays
disp("    Cell Arrays");
c = {[10 20 30],[40],[50 60]}
c{1}
c{1}(1)
c(end+1) = {1:2:10}
c{end}
t = {'Miami','Detroit','Boston'}
t = sort(t)
x = [5.0895 1.9663];
P = [1 1; 6 1; 6 5];
xP = {x, P}
d = mydist(xP{:})
c = cell(1,3)
[c2{1:3}] = deal(0)
[c3{1:3}] = deal(1,2,3)
%% Structures
disp("    Structures");
s.Name = 'Mike'
s.Age = 44
s(2).Name = 'Bill'
s(2).Age = 40;
s(2)
s.Name
s = struct('Name',{'Mike','Bill'},'Age',{44,40})
names = {s.Name}
ages = [s.Age]
s2.Name = {'Mike','Bill'}
s2.Age = [44 40]
%% N-D Arrays
disp("    N-D Arrays");
D3 = [1 2 3; 4 5 6]
D3(:,:,2) = [7 8 9; 10 11 12]
D3(1,3,2)
D3(:,1,:)
D2 = squeeze(D3(:,1,:))
%% 13. Control Structures
disp(" 13. Control Structures");
%% FOR Loop
disp("    FOR Loop");
for i = 1:3
i
end
for i = 5:-2:1, i, end
chararray = 'abc'
for i = chararray
i
end
c = {[10 20 30],[40],[50 60]}
c = c + 1
for i = 1:length(c), c{i} = c{i} + 1;
end
c{:}
%% IF Statement
disp("    IF Statement");
n = 3;
if n > 0
disp('Positive value.')
elseif n < 0
disp('Negative value.')
else
disp('Zero.')
end
%% WHILE Loop
disp("    WHILE Loop");
while n > 0
n = n - 1
end
%% DO-WHILE Loop
disp("    DO-WHILE Loop");
done = false;
while ~done
n = n + 1
if n >= 3, done = true; end
end
a = [5 0 -1 9 0];
a > 0
any(a > 0)
all(a > 0)
%% 14. Example: Random Walk Simulation
disp(" 14. Example: Random Walk Simulation");
rng(123)
s = rand(1,5)
s > 0.5
(s > 0.5)*2
((s > 0.5)*2) - 1
d = cumsum(((s > 0.5)*2) - 1)
s = rand(1,100);
d = cumsum(((s > 0.5)*2) - 1);
plot(d)
for i = 1:100
d = cumsum(((rand(1,1000)>0.5)*2)-1);
dt(i) = d(end);
end
mean(abs(dt))
%% 15. Logical vs. Index Arrays
disp(" 15. Logical vs. Index Arrays");
a
ispos = a > 0
a(ispos)
idxpos = find(a > 0)
a(idxpos)
s = {'Miami','Detroit','Boston'};
idxDetroit = strmatch('Detroit',s)
isDetroit = strcmpi('detroit',s)
%% Advantages of Index Arrays
disp("    Advantages of Index Arrays");
a
[sa, idxa] = sort(a)
idx = [1 2 1];
a(idx)
[mina, idxmina] = min(a)
%% 2. Referring to all of cities in North Carolina and South Carolina
% with populations of at least 10,000:
disp('Question 2');
%% (a) Assuming weight is proportional to population, what location minimizes
% the sum of the weight times the great-circle distance to each city?
s = uscity(strcmp('NC',uscity('ST')) & uscity('Pop') >= 10000 | strcmp('SC',uscity('ST')) & uscity('Pop') >= 10000);
TCh = @(xy) dists(xy,s.XY,'mi') * s.Pop;
xy = fminsearch(TCh,mean(s.XY));
fprintf("(a) Minisum location for NC and SC is %f, %f.\n" , xy(1), xy(2));
%% (b) What city is closest to the optimal location?
x = s.Name{argmin(dists(xy,s.XY,'mi'))};
fprintf("(b) Closest city is %s, %s.\n" , x, s.ST{strcmp(s.Name, x)})
%% (c) What city is furthest from the location?
[n, idx] = sort(dists(xy,s.XY,'mi'));
x = strjoin(s.Name(idx(end)));
fprintf("(c) The furthest city from the location is %s, %s.\n" , ...
    char(x), char(s.ST(strcmp(s.Name, x))))
%% (d) What are the four closest cities to the location?
x = {s.Name{idx(1:4)}};
x = strjoin(x, ", ");
fprintf("(d) The four closest cities to the location are %s.\n" , x)
%% (e) Plot the name and location of the four closest cities.
x = s.XY(idx(1:4),:);
makemap([x; xy]);
xyc = {'Optimal Point'};
pplot(xy, 'g*');
pplot(xy, xyc);
pplot(x,'r.');
pplot(x,s.Name(idx(1:4)));
%% (f) How far is the location from the largest (greatest population) city?
bgst = argmax(s.Pop);
x = dists(xy,s.XY(bgst,:),'mi');
name = s.Name{bgst};
fprintf("(f) The location is %s miles away from %s.\n" , num2str(x), name);
%% (g) What city with a population of at least 50,000 is closest to the
% location?
m = s.Pop >= 50000;
idxpos = find(m > 0);
x = s.Name{idxpos(argmin(dists(xy,s.XY(idxpos,:),'mi')))};
fprintf("(g) The closest city with a population of at least 50,000 is %s.\n" , x)
%% (h) What percentage of the population in the cities is South of the
% location?
idxpos = find(s.XY(:,2)< xy(2));
x = sum(s.Pop(idxpos))/sum(s.Pop)*100;
fprintf("(h) Percentage of the population in the cities South of the location is %.3f%%.\n" , x)
%% (i) What is the total population of all of the cities that are within
% 100 miles of the location?
idxpos = find(dists(xy,s.XY,'mi') < 100);
x = sum(s.Pop(idxpos));
fprintf("(i) The total population of all of the cities that are within 100 miles of the location is %d people.\n" , x)
%% (j) Plot the name and location of the five closest cities to the
% optimal location that have a population of at least 50,000.
idxpos = find(m > 0);
x = s.XY(idxpos,:);
name = s.Name(idxpos);
[n idx] = sort(dists(xy,x,'mi'));
x = x(idx,:);
name = name(idx);
makemap([x(1:5,:); xy]);
pplot(xy, 'g*');
pplot(xy,xyc);
pplot(x(1:5,:),'r.');
pplot(x(1:5,:),name(1:5));
%% (k) Assuming that three DCs are located in Raleigh, NC, Charlotte, NC,
% and North Charleston, SC, respectively, and that each city is served by
% its closest DC, what is the total population served by each DC?
city2lonlat = @(city,st) uscity('XY',mand(city,uscity('Name'),st,uscity('ST')));
P = city2lonlat({'Raleigh','Charlotte','North Charleston'},{'NC','NC','SC'});

hold on
voronoi(P(:,1),P(:,2));
pplot(s.XY,'r.');
pplot(s.XY,s.Name), shg;

[vx, vy] = voronoi(P(:,1),P(:,2));
mh = @(i)(vy(1)-vy(i))/(vx(1)-vx(i)); % slope of the line
bh = @(i) vy(1)- mh(i)*vx(1);         % y-intercept
y0 = @(i,x) (mh(i)*x(1)+bh(i));       % line equation

NorChar = [];
Ral = [];
Char = [];
for i = 1:length(s.Name)
    if s.XY(i,1) > vx(1) % right to the central point
        if s.XY(i,2) < y0(6,s.XY(i,1)) % closer to North Charleston
        NorChar(end+1) = i;
        else
            Ral(end+1) = i;
        end
    elseif s.XY(i,2) < y0(4,s.XY(i,1)) % left to the central point and
        NorChar(end+1) = i;            % closer to North Charleston
    elseif s.XY(i,2) > y0(2,s.XY(i,1)) % closer to Raleigh
        Ral(end+1) = i;
        else Char(end+1) = i;
    end
end
Raleigh = sum(s.Pop(Ral));
Charlotte = sum(s.Pop(Char));
NCharleston = sum(s.Pop(NorChar));
sum(s.Pop)- sum([Raleigh Charlotte NCharleston]); %checking if the answer is correct
fprintf("(k) the total population served by each DC");
vdisp('Raleigh, Charlotte, NCharleston');
%% 3. Assume that 20, 30, and 24 full truckloads (TL) per year of finished
% goods will be shipped to customers located in Detroit, MI, Gainesville, FL,
% and Memphis, TN.
disp('Question 3')
city2lonlat = @(city,st) uscity('XY',mand(city,uscity('Name'),st,uscity('ST')));
P = city2lonlat({'Detroit','Gainesville','Memphis'},{'MI','FL','TN'});
f = [20 30 24];
%% (a) Where (in decimal degrees) should a new factory be located in order
% to minimize total outbound truck travel?
xy = minisumloc(P,f,'mi');
fprintf("(a) A new factory should be located at %f, %f.\n", xy(1), xy(2));
%% (b) Near what city with a population of at least 50,000 is this location?
s = uscity50k;
x = [s.Name{argmin(dists(xy,s.XY,'mi'))} s.ST(strcmp(s.Name,...
    s.Name{argmin(dists(xy,s.XY,'mi'))}))];
fprintf("(b) The closest city with a population of at least 50,000 is %s.\n"...
    , strjoin(x, ", "))
%% (c) Near what city with a population of at least 10,000 is this location?
s = uscity10k;
x = [s.Name{argmin(dists(xy,s.XY,'mi'))} s.ST(strcmp(s.Name,...
    s.Name{argmin(dists(xy,s.XY,'mi'))}))];
fprintf("(c) The closest city with a population of at least 10,000 is %s.\n"...
    , strjoin(x, ", "));
%% (d) A circuity factor is the average ratio of the actual road distance between
% two points and the great circle distance between the two points. Use the 
% simple average of the distances between all pairs of customers (three
% values, in total) to estimate the circuity factor. You can use Google Map
% to determine the actual road distances (where you can pick the fastest 
% route if several routes are available).
D = dists(P,P,'mi');
g = mean([1053.7/D(1,2) 744/D(1,3) 721.8/D(2,3)]);
fprintf("(d) A circuity factor is %.3f.\n", g);
%% (e) If the transportation rate is $3.00 per mile for each one-way TL
% shipment, what is the total transportation cost per year assuming that
% the trucks do not return to the factory and that the actual road 
% distances are estimated using a circuity factor?
r = 3;
[xy,TC] = minisumloc(P,f*r*g,'mi');
fprintf("(e) The total transportation cost per year is %.2f $.\n", TC);
%% (f) Assuming that the factory will be built in Cary, NC instead of at 
% the optimal location, what is the increase in total transportation cost
% per year?
xy = city2lonlat('Cary','NC');
TC_Cary = sum(f.*dists(xy,P,'mi')) * g * r;
vdisp('TC,TC_Cary,TC_Cary - TC')
fprintf("(f) The increase in total transportation cost per year is %.2f $.\n",...
    TC_Cary - TC);
