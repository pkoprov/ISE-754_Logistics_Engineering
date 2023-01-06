function TC = myrteTC(rte,rTDh0,maxtime)
% My VRP constraints
TC = rTDh0(rte);
if TC > maxtime
   TC = Inf;
end
