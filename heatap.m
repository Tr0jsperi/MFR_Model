function [cpa] = heatap(tp)
% c calculates ash heat capacity
% c tp particle temperature (k)
% c cpa ash heat capacity (cal/g/k)
% c rgas gas constant (1.987 cal/g/k)
cpa = (754+0.586*(tp-273))/(4.184e3);