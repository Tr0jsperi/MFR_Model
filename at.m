function [tri_area] = at(aa,bb,cc)
% Function required by lightgas.m to calculate area of triangle
tri_area = 0.5*bb*cc*(1-((bb^2+cc^2-aa^2)/(2*bb*cc))^2)^0.5;