function r = myroots_cmm(a)
% MYROOTS Find the roots of a polynomial using the companion matrix method
%   R = MYROOTS(A) finds the roots of the polynomial with coefficients A using
%   the companion matrix method. The input A is a vector of length N+1 containing
%   the coefficients of the polynomial in descending order:
%
%     p(x) = a(1)*x^N + a(2)*x^(N-1) + ... + a(N)*x + a(N+1)
%
%   The output R is a vector of length N containing the roots of the polynomial.

n = length(a) - 1; % Degree of polynomial
C = zeros(n);
C(1,:) = -a(2:end) / a(1);
C(2:end,1:n-1) = eye(n-1);
r = eig(C);
end

