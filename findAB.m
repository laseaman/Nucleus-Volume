function [ alpha, beta, F, Fs ] = findAB( fi,w,times )
%Finds Alpha and Beta for fitting a single frequency pair of basis
%functions.
%   Used by FreqBasis.m, Finds a,b to minimize F=sum(a*si+b*ci-fi)^2. fi is
%   a vector of the data, w is the frequency to be fit to it. 
    si = sin(2*pi*w.*times); ci = cos(2*pi*w.*times);
    A = [nansum(si.^2) nansum(ci.*si); nansum(ci.*si) nansum(ci.^2)];
    a_top = [nansum(fi.*si) nansum(ci.*si); nansum(fi.*ci) nansum(ci.^2)];
    b_top = [nansum(si.^2) nansum(fi.*si); nansum(ci.*si) nansum(fi.*ci)];
    alpha = det(a_top)/det(A);
    beta = det(b_top)/det(A);
    F = nansum((alpha.*si+beta.*ci-fi).^2);
    Fs = nansum(fi.^2) - nansum( (alpha.*si+beta.*ci-fi).^2);
end
