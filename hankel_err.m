function [err,tmp] = hankel_err(G,Gp,tol)
    err_curr = -Inf; rate = 0.02;
    zs = sqrt(3):rate:(2*pi+sqrt(3));
    zs = [zs,0:rate:2*pi];
    zs = cos(zs) + 1j*sin(zs);
    tmp = [];
    for i = 1:length(zs)
        err = norm(G(zs(i))-Gp(zs(i)));
        err_curr = max([err_curr,err]);
        if err > 0.004
            tmp = [tmp,i];
        end
    end
    err = err_curr;
end
