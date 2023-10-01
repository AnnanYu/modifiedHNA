% Modified HNA algorithm
function [Ahat,Bhat,Chat,Dhat,t,t1] = modified_HNA(A,B,C,D,r,tol,Utol)
    rng(37);
    tic
    n = length(A); [p,m] = size(D);
    if r >= n
        Ahat = A; Bhat = B; Chat = C; Dhat = D; t = toc; t1 = 0;
        return
    end
    sys = ss(A,B,C,D);
    [sysb,sigma] = balreal(sys);
    A = sysb.A; B = sysb.B; C = sysb.C; D = sysb.D;
    sighat = sigma(r+1); ind2 = find(abs(sigma - sighat) < tol);
    ind1 = 1:(ind2(1)-1); ind3 = (ind2(length(ind2))+1):n;
    A = [A(ind1,ind1),A(ind1,ind3),A(ind1,ind2);...
        A(ind3,ind1),A(ind3,ind3),A(ind3,ind2);...
        A(ind2,ind1),A(ind2,ind3),A(ind2,ind2)];
    B = [B(ind1,:);B(ind3,:);B(ind2,:)];
    C = [C(:,ind1),C(:,ind3),C(:,ind2)];
    sigma = [sigma(ind1);sigma(ind3);sigma(ind2)]; P = diag(sigma);
    t = toc;
    
    cut = n - length(ind2);
    ind1 = 1:cut; ind2 = (cut+1):n;
    Sighat = diag((sigma(ind1).^2-sighat^2).^-1);
    X = -C(:,ind2)'; Y = B(ind2,:);
    [U,S,V] = svd(X); rr = length(find(diag(S) > Utol));
    Vnew = (U(:,1:rr)*S(1:rr,1:rr)) \ Y; Vnew = Vnew';
    for j = 1:(p-rr)
        vrand = rand(m,1);
        vortho = vrand - Vnew * (Vnew' * vrand);
        vnormal = vortho / norm(vortho);
        Vnew = [Vnew,vnormal];
    end
    U = V * Vnew';
    % U = -C2 * pinv(B2');

    Ahat = Sighat * (sighat^2 * A(ind1,ind1)' ...
        + P(ind1,ind1)*A(ind1,ind1)*P(ind1,ind1) ...
        - sighat * C(:,ind1)' * U * B(ind1,:)');
    Bhat = Sighat * (P(ind1,ind1) * B(ind1,:) + sighat * C(:,ind1)' * U);
    Chat = C(:,ind1) * P(ind1,ind1) + sighat * U * B(ind1,:)';
    Dhat = D - sighat * U;
    t1 = toc; t1 = t1 - t;
end