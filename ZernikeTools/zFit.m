function [zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,z,maxDegree)
% [zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,maxDegree)
%
% Least-squares fit of Zernike polynomials to data.
%
% rho = column vector of normalized radius
% theta = column vector of angles (rad)
% z = col vector of data
% maxDegree = maximum degree of radial polynomials. All polynomials through 
% maxDegree are fit.
%
% zCoef = col vector of fitted coefficients (same units as z)
% zMatrix = matrix whose columns are the Zernike polynomials evaluated over
% the data (useful for subtracting the fitted polynomials)
% nVec = col vector of radial degrees of polynomials
% elVec = col vector of angular orders of polynomials
%
% Polynomials are normalized to +/-1 at edge of pupil.
%
[zMatrix,nVec,elVec] = zEval(rho,theta,maxDegree);
zCoef = zMatrix \ z;
