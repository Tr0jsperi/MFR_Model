function [dist] = d(x2,x1,y2,y1)
% Function required by lightgas.m to calculate distance
dist = ((x2-x1)^2+(y2-y1)^2)^0.5;