function [A,B,C,t,F] = aaa_reduction(G,m,p,samples,d,tol,svd_tol)
    Z = cos(samples)'+1j*sin(samples)';
    F = zeros(length(Z),m*p);
    for i = 1:length(Z)
        F(i,:) = reshape(G(Z(i)),1,m*p);
    end
    
    tic
    [~,pol,~,~,~,errvec] = aaa2(F,Z,tol,d);
    deg = length(pol); coef = zeros(p,m,deg);
    M = zeros(length(Z),deg);
    for i = 1:deg
        M(:,i) = (Z-pol(i)).^(-2);
    end
    C = M \ F;
    for i = 1:p
        for j = 1:m
            coef(i,j,:) = C(:,(i-1)*m+j);
        end
    end
    A = []; B = []; C = [];
    AF = []; BF = []; CF = [];
    for i = 1:deg
        if abs(pol(i)) > 1 - 1e-8
            slice = transpose(coef(:,:,i)); [U,S,V] = svd(slice);
            s = diag(S); s = s(s > svd_tol); r = length(s);
            AF = blkdiag(AF,pol(i)*eye(r)); BF = [BF,V(:,1:r)*S(1:r,1:r)];
            CF = [CF,U(:,1:r)];
            continue
        end
        slice = transpose(coef(:,:,i)); [U,S,V] = svd(slice);
        s = diag(S); s = s(s > svd_tol); r = length(s);
        A = blkdiag(A,pol(i)*eye(r)); B = [B,V(:,1:r)*S(1:r,1:r)];
        C = [C,U(:,1:r)];
    end
    B = B';
    t = toc;
    F = @(z) CF * ((z*eye(length(AF)) - AF) \ BF');
    if isempty(F(1))
        F = @(z) 0;
    end
    disp(length(A));
end