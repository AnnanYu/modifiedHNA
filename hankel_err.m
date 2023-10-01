function [err,tmp] = hankel_err(G,Gp,tol)
    %err_prev = 0;
    err_curr = -Inf; rate = 0.02;
    %zs = [logspace(0.5,5.5,300)*1j,-logspace(0.5,5.5,300)*1j];
    %zs(abs(zs - pi) < 0.005) = [];
    %zs = [2.88:0.000005:2.90,3.04:0.000005:3.06,...
    %    3.22:0.000005:3.24,3.38:0.000005:3.40];
    zs = sqrt(3):rate:(2*pi+sqrt(3));
    %zs = 0:rate:(2*pi);
    zs = [zs,0:rate:2*pi];
    zs = cos(zs) + 1j*sin(zs);
    %while (abs(err_prev - err_curr) > tol)
    %    err_prev = err_curr;
    tmp = [];
    for i = 1:length(zs)
        err = norm(G(zs(i))-Gp(zs(i)));
        err_curr = max([err_curr,err]);
        if err > 0.004
            tmp = [tmp,i];
        end
    end
    %    rate = rate / 3;
    %    zs = 0:rate:2*pi; zs = cos(zs) + 1j*sin(zs);
    %end
    err = err_curr;
end