% Test the accuracy of the modified HNA algorithm

%% Setup
rng(37);
d = rand(16,1);
X = rand(16,16);
X = orth(X);
X = X * diag(d) * X';

% Construct a grid
delta = 10.^(-1.1:-0.5:-11); theta = 10.^(-1:-0.5:-11);
errs = zeros(length(delta),length(theta));

%% For each configuration of delta and theta, run the modified HNA algorithm
for i = 1:length(delta)
    for j = 1:length(theta)
        [A,B,C,P] = system_factory(X,theta(j),delta(i));
        [Ahat,Bhat,Chat,Dhat] = glover_optimal(A,B,C,diag(P));
        [A,B,C,D] = plane2disk(A,B,C,zeros(16,16));
        [Ahat,Bhat,Chat,Dhat] = plane2disk(Ahat,Bhat,Chat,Dhat);
        G = @(z) D + C * ((z*eye(length(A)) - A) \ B);
        Ghat = @(z) Dhat + Chat * ((z*eye(length(Ahat)) - Ahat) \ Bhat);
        errs(i,j) = hankel_err(G,Ghat,theta) - 0.5;
    end
end

%% Plot the results
figure;
X = delta' * ones(1,21); Y = ones(20,1)*theta;
surf(X,Y,errs)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ZScale','log')
set(gca,'ColorScale','log')
colorbar

figure;
mat = abs(Ahat);
X = ones(11,1)*(1:11);
Y = X';
surf(X,Y,mat)
set(gca,'ZScale','log')
set(gca,'ColorScale','log')
colorbar

%% Modified HNA algorithm
function [Ahat,Bhat,Chat,Dhat] = glover_optimal(A,B,C,sigma)
    ind1 = 1:11; ind2 = 12:16; sighat = 0.5; n = 16; P = diag(sigma);
    Sighat = diag((sigma(ind1).^2-sighat^2).^-1);
    X = -C(:,ind2)'; Y = B(ind2,:);
    [U,S,V] = svd(X); rr = length(find(diag(S) > 1e-3));
    Vnew = (U(:,1:rr)*S(1:rr,1:rr)) \ Y; Vnew = Vnew';
    for j = 1:(n-rr)
        vrand = rand(n,1);
        vortho = vrand - Vnew * (Vnew' * vrand);
        vnormal = vortho / norm(vortho);
        Vnew = [Vnew,vnormal];
    end
    U = V * Vnew';

    Ahat = Sighat * (sighat^2 * A(ind1,ind1)' ...
        + P(ind1,ind1)*A(ind1,ind1)*P(ind1,ind1) ...
        - sighat * C(:,ind1)' * U * B(ind1,:)');
    Bhat = Sighat * (P(ind1,ind1) * B(ind1,:) + sighat * C(:,ind1)' * U);
    Chat = C(:,ind1) * P(ind1,ind1) + sighat * U * B(ind1,:)';
    Dhat = -sighat * U;
end