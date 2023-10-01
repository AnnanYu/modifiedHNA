% Approximate the Hilbert operator using the AAA algroithm
% and the modified HNA algorithm.

%% Create a transfer function
G = @(z) sum((1:500).^(-1) .* ((z).^(-1:-1:-500)));

%% Sample the transfer function
deg = 10:5:70; es = zeros(size(deg)); ts = zeros(3,length(deg));
theta = [logspace(-4,0,100)-(1e-4),(1e-4)-1*logspace(-4,0,100),0:0.01:(2*pi)];
theta(theta == 0) = []; theta = [theta,0];

samples = cos(theta) + 1j*sin(theta);
FF = cell(length(samples),1);
for i = 1:length(samples)
    FF{i} = G(samples(i));
end

%% Run block-AAA and modified HNA
opts.reg = 0;
for i = 1:length(deg)
    opts.maxit = deg(i);
    [R1,~,~,Ap,Bp,Cp,t1,Hinter,Ginter,maxerr] = block_aaa_disk(FF,samples,opts);
    [Ap,Bp,Cp,Dp] = disk2plane(Ap,Bp,Cp,0);
    [Ap,Bp,Cp,Dp,t2,t3] = modified_HNA(Ap,Bp,Cp,Dp,10,0.00001,1e-10);
    [Ap,Bp,Cp,Dp] = plane2disk(Ap,Bp,Cp,Dp);
    Gp = @(z) Dp + Cp * ((z*eye(length(Ap)) - Ap) \ Bp);
    es(i) = hankel_err(G,@(z)Gp(z)+Hinter(z),1e-4);
    ts(1,i) = t1; ts(2,i) = t2; ts(3,i) = t3;
end

%% Plot the results
figure; hold on
plot(deg,es,'-x','linewidth',2,'markersize',8,'displayname','error');
plot([deg(1),deg(length(deg))],[0.0360,0.0360],'--','linewidth',3,...
    'displayname','theo. low. bd.');
xlabel('degree of AAA')
ylabel('Hankel error')
legend show
grid on
hold off

figure; hold on
plot(deg,ts(1,:),'-x','linewidth',2,'markersize',8,'displayname','t1');
plot(deg,ts(2,:),'-x','linewidth',2,'markersize',8,'displayname','t2');
plot(deg,ts(3,:),'-x','linewidth',2,'markersize',8,'displayname','t3');
plot(deg,ts(1,:)+ts(2,:)+ts(3,:),'-x','linewidth',2,...
    'markersize',8,'displayname','total time');
xlabel('degree of AAA')
ylabel('time')
legend show
grid on
hold off
