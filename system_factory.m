function [A,B,C,P] = system_factory(X,theta,delta)
sigma = [0.1:0.1:0.4,0.6:0.1:1];
% sigma = [sigma,0.5,0.5,0.5,0.5,0.5,0.5]; sigma = sort(sigma); sigma = [sigma,0.5+delta];
sigma = [sigma,0.5-delta,0.5+delta];
sigma = sort(sigma);
sigma = [sigma,0.5-theta,0.5-theta/2,0.5,0.5+theta/2,0.5+theta];

P = diag(sigma); A = -X;
B = chol(-A*P-P*A'); B = B'; C = chol(-A'*P-P*A);
end