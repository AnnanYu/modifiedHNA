% Test the advantage of the modified HNA algorithm
% on the iss example

%% Load the matrices
clear all
rng(37);
load('benchmarks/iss.mat');
D = zeros(3,3);
[Ar,Br,Cr,Dr,~,~] = modified_HNA(A,B,C,D,175,1e-15,0);
[Amr,Bmr,Cmr,Dmr,~,~] = modified_HNA(A,B,C,D,175,1e-12,0);

%% Test it with a single epsilon
[Ar,Br,Cr,Dr] = plane2disk(Ar,Br,Cr,Dr);
[Amr,Bmr,Cmr,Dmr] = plane2disk(Amr,Bmr,Cmr,Dmr);
[Ad,Bd,Cd,Dd] = plane2disk(A,B,C,D);
iter = 10000; u = rand(3,iter);
y = simulate(Ad,Bd,Cd,Dd,u,iter);
yr = simulate(Ar,Br,Cr,Dr,u,iter);
ymr = simulate(Amr,Bmr,Cmr,Dmr,u,iter);
fprintf('Error of HNA: %.3e\n',norm(y-yr,'fro')/sqrt(iter));
fprintf('Error of modified HNA: %.3e\n',norm(y-ymr,'fro')/sqrt(iter));

%% Test it with multiple epsilons
dists = abs(hsv - hsv(175)); dists = sort(dists);
es = (dists(2:length(dists)) + dists(1:(length(dists) - 1))) / 2;
es = sort(es); es(es > 1e-2) = [];
errs = zeros(size(es)); maxent = zeros(size(es));
for i = 1:length(es)
    [Amr,Bmr,Cmr,Dmr,~,~] = modified_HNA(A,B,C,D,175,es(i),0.01);
    maxent(i) = max(max(abs(Amr)));
    [Amr,Bmr,Cmr,Dmr] = plane2disk(Amr,Bmr,Cmr,Dmr);
    errs(i) = norm(y-simulate(Amr,Bmr,Cmr,Dmr,u,iter),'fro')/sqrt(iter);
    if i == 120
        disp(i)
    end
end

%% Plot the results
es = zeros(2*(length(errs)+1));
for i = 1:length(errs)
    es(2*i-1) = dists(i); es(2*i) = dists(i);
end
es = es(2:(length(es)-1))'; es(1) = 10^-13;
err = zeros(2*length(errs),1);
for i = 1:length(errs)
err(2*i-1) = errs(i); err(2*i) = errs(i);
end
errs = err;
es(isnan(errs)) = []; maxent(isnan(errs)) = []; errs(isnan(errs)) = [];
figure;
loglog(es,errs,'linewidth',2);
hold on
plot([0.008,0.00001],[0.0001,1.25e-7],'--','linewidth',1.5)
grid on
box on
axis([10^-13,0.006,10^-10,10^10])


%% Helper function, simulate the system
function y = simulate(A,B,C,D,u,iter)
    rng(37);
    x = zeros(length(A),1);
    y = zeros(3,iter);
    for i = 1:iter
        x = A * x + B * u(:,i);
        y(:,i) = C * x + D * u(:,i);
    end
end