function [zMatrix,nVec,elVec] = zEval(rho,theta,maxDegree)
% [zMatrix,nVec,elVec] = zEval(rho,theta,maxDegree)
%
% Evaluates Zernike polynomials over the points specified by col vectors
% rho and theta.
%
% rho = column vector of normalized radius
% theta = column vector of angles (rad)
% maxDegree = maximum degree of radial polynomials. All polynomials through 
% maxDegree are evaluated.
%
% zMatrix contains polynomials as its columns.
% nVec = col vector of radial degrees of polynomials
% elVec = col vector of angular orders of polynomials
%
% Polynomials are normalized to +/-1 at edge of pupil.
%
% Use notation of Born & Wolf, Appendix VII, with the additional convention
% that el>0 represents cosine term and el<0 represents sine term. (el
% because I can't tell an l from a 1)
% Put polynomials in the order in which they appear in Durango, which uses
% m for what I call el.
%
points = length(rho);
i = 1;  % index for polynomials
for n = 0:maxDegree % degree of radial polynomial
    for el = n:-2:-n
        m = abs(el); % angular order
        radialPoly = zeros(points,1);
        for s = 0:(n-m)/2
            coef = (-1)^s*factorial(n-s) / ...
                (factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s));
            power = n-2*s;
            radialPoly = radialPoly + coef * rho.^power;
        end
        if el>0
            zMatrix(:,i) = radialPoly .* cos(m*theta);
            nVec(i,1) = n;
            elVec(i,1) = el;
        elseif el==0
            zMatrix(:,i) = radialPoly;
            nVec(i,1) = n;
            elVec(i,1) = el;
        else
            zMatrix(:,i) = radialPoly .* sin(m*theta);
            nVec(i,1) = n;
            elVec(i,1) = el;
        end
        i = i+1;
    end
end
