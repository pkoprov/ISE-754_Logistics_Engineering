%% HW1
%% 1
disp("Problem 1")
x = [3 1 2 9 5 4]
%% (a) Extract the third element from x
disp("(a) Extract the third element from x")
x(3)
%% (b)Extract all but the last element from x
disp("(b)Extract all but the last element from x")
x(1:end-1)
%% (c) Extract the first, third, first, sixth, and first element from x
disp("(c) Extract the first, third, first, sixth, and first element from x")
x([1 3 1 6 1])
%% (d) Reverse the elements of x
disp("(d) Reverse the elements of x")
fliplr(x)
%% (e) Calculate the sum of all of the elements of x
disp("(e) Calculate the sum of all of the elements of x")
sum(x)
%% (f) Calculate the sum from the first to the ith element of x, for all 
%      elements 1 to i in x
disp(sprintf("(f) Calculate the sum from the first to the ith element of x, for all" + "\n" + "elements 1 to i in x"))
cumsum(x)
%% (g) Modify x by setting the second and sixth elements of x equal to zero
disp("(g) Modify x by setting the second and sixth elements of x equal to zero")
x([2 6]) = 0
%% (h) Using the result from (g), modify x by deleting its third element
disp("(h) Using the result from (g), modify x by deleting its third element")
x(3) = []
%% (i) Using the result from (h), modify x by adding 7 to its end
disp("(i) Using the result from (h), modify x by adding 7 to its end")
x = [x 7]
%% (j) Using the result from (i), modify x by converting it into a 2 x 3
%      matrix, where the first row of the matrix has the first three
%      elements of x
disp(sprintf("(j) Using the result from (i), modify x by converting it \n into a 2 x 3 matrix, where the first row of the matrix has \n the first three elements of x"))
x = [x(1:3); x(4:6)]
%% 2
disp ("Problem 2")
x = [6 2 1 4]
A = [2 7 9 7; 3 2 5 6; 8 2 1 5]
%% (a) Add x to each row of A
disp("(a) Add x to each row of A")
A + x
%% (b) Add x to the sum of each column of A
disp("(b) Add x to the sum of each column of A")
sum(A)+ x
%% (c) Add twice the sum of x to each element of A
disp("(c) Add twice the sum of x to each element of A")
2*sum(x) + A
%% (d) Calculate the element-by-element product of each row of A and x
disp("(d) Calculate the element-by-element product of each row of A and x")
A .* x
%% (e) For each row of A, calculate the sum of the element-by-element
%      product of each row and x
disp(sprintf("(e) For each row of A, calculate the sum of the element-by-element\n product of each row and x"))
sum(A .*x,2)
%% 3
disp("Problem 3")
x = [1 5 2 7 9 0 1]
y = [5 1 2 8 0 0 2]
%% (a) Extract from x those values that are greater than the
%      corresponding values of y
disp(sprintf("(a) Extract from x those values that are greater than the\n corresponding values of y"))
x(x > y)
%% (b) Extract from x those values that are both greater than the
%      corresponding values of y and less than 6
disp(sprintf("(b) Extract from x those values that are both greater than the \n corresponding values of y and less than 6"))
x((x > y) & (x < 6))
%% (c) Extract from x those values that are either less than 2 or 
%      greater than 6
disp(sprintf("(c) Extract from x those values that are either less than 2 or \n greater than 6"))
x((x < 2) | (x > 6))
%% (d) Modify y by adding 1 to each of its nonzero values
disp("(d) Modify y by adding 1 to each of its nonzero values")
y(y > 0) = y(y > 0) + 1
%% (e) Divide each element of y by the corresponding element of x 
%      as long as the element of x is nonzero (to avoid dividing by zero)
disp(sprintf("(e) Divide each element of y by the corresponding element of x \n as long as the element of x is nonzero (to avoid dividing by zero)"))
y(x~=0)./x(x~=0)
%% (f) Modify y by setting all of its zero values to 1
disp("(f) Modify y by setting all of its zero values to 1")
y(y == 0) = 1
%% 4
disp("Problem 4")
% Modify the Minimum-Distance Location example in Basic Concepts so that
% it can be used to find the location that minimizes the maximum distance
% traveled between x and the three points in P. (Note: you need to create
% the function mydist discussed in the example.)
disp(sprintf("Modify the Minimum-Distance Location example in Basic Concepts so that\nit can be used to find the location that minimizes the maximum distance\ntraveled between x and the three points in P. (Note: you need to create\nthe function mydist discussed in the example.)"))
x = [3 1];
P = [1 1; 6 1; 6 5];
MaxMydist = @(x) max(mydist(x,P));
[x,maxd] = fminsearch(MaxMydist,[0 0])

%% 5
disp ("Problem 5")
% Modify the Minimum-Distance Location example in Basic Concepts so that
% it can be used to find the location that minimizes the sum of distance
% traveled assuming that 3, 4, and 2 trips are made between x and the
% three points in P, respectively.
disp (sprintf("Modify the Minimum-Distance Location example in Basic Concepts so that\nit can be used to find the location that minimizes the maximum distance\ntraveled between x and the three points in P. (Note: you need to create\nthe function mydist discussed in the example.)"))
x = [3 1];
P = [1 1; 6 1; 6 5];
t = [3 4 2]'
SumMydist = @(x) sum(t .* mydist(x,P));
[x,sumd] = fminsearch(SumMydist, x)
%% 6
disp("Problem 6")
% Write a code to generate the first 10 numbers of the Fibonacci Sequence.
disp("Write a code to generate the first 10 numbers of the Fibonacci Sequence.")
n = 1:10;
disp(sprintf("\nBuilt-in function 'fibonacci' "))
fibonacci(n)
% or
disp("Or own script")
x = 1;
y = [1 1];
for i = 1:8
    x = sum(y(end-1 : end));
    y = [y x];
    if i == 8
        y
    end
end
%% 7
disp("Problem 7")
% Plot the quadratic function, 𝑓(𝑥)=𝑥2+2𝑥+1, for 𝑥∈[−20,20].
disp("Plot the quadratic function, 𝑓(𝑥)=𝑥2+2𝑥+1, for 𝑥∈[−20,20].")
x = -20:20;
y = x.^2 + 2*x + 1;
plot(x, y)
%% 