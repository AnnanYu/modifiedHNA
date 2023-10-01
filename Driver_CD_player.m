% Driver for reducing the CD_Player example

%% Load the matrices
clear all
index = 2;

samples1 = 0:0.005:2*pi; samples1 = cos(samples1) + 1j * sin(samples1);
samples1 = (samples1-1) ./ (samples1+1);
samples3 = 0:0.001:2*pi;
samples3 = cos(samples3) + 1j * sin(samples3);
samples3 = (samples3-1) ./ (samples3+1);
samples2 = [logspace(1,3,20)*1j,-logspace(1,3,20)*1j];
configs = {{'benchmarks\eady.mat',10,10:1:40,samples1,0.01,1},...
    {'benchmarks\CDplayer.mat',10,10:2:60,samples2,0.01,1},...
    {'benchmarks\iss.mat',10,20:5:80,samples3,0.01,1}};
config = configs{index};
load(config{1})

% Compute missing variables
if ~exist('C','var')
    C = B';
end
[p,~] = size(C); [~,m] = size(B);
if ~exist('D','var')
    D = zeros(p,m);
end
if ~exist('hsv','var')
    sys = ss(A,B,C,D);
    [~,hsv] = balreal(sys);
end

% Get the target degree, AAA degrees, and benchmarks
target_d = config{2};
err_bench = hsv(target_d+1); % Theoretical lower bound
deg = config{3};
[~,~,~,~,t2,t3] = modified_HNA(A,B,C,D,target_d,1e-10,1e-10);
t_bench = t2 + t3;
samples = config{4};
err_conf = config{5};
demo = config{6};

G1 = @(z) C * ((z*eye(length(A)) - A) \ B);
FF = cell(length(samples),1);
for i = 1:length(samples)
    FF{i} = G1(samples(i));
end

%% First stage of AAA
[Ap,Bp,Cp,Dp] = plane2disk(A,B,C,D);
G = @(z) Dp + Cp * ((z*eye(length(Ap)) - Ap) \ Bp);
sys = ss(A,B,C,D);
[sysb,hsv] = balreal(sys);
Ab = sysb.A; Bb = sysb.B; Cb = sysb.C; Db = sysb.D;
Ab = Ab(1:target_d,1:target_d); Bb = Bb(1:target_d,:); Cb = Cb(:,1:target_d);
[Ab,Bb,Cb,Db] = plane2disk(Ab,Bb,Cb,Db);
Gb = @(z) Db + Cb * ((z*eye(length(Ab)) - Ab) \ Bb);
errb = hankel_err(G,Gb,err_conf); % Balanced truncation error bound

opts.maxit = 10; opts.reg = 1e-2;
[R1,~,~,A1,B1,C1,t1,Hinter,Ginter,maxerr] = block_aaa(FF,samples,opts);
samples = [logspace(0.01,4,1000)*1j,-logspace(0.01,4,1000)*1j];
for i = 1:length(samples)
    FF{i} = G1(samples(i)) - R1(samples(i));
end

%% Second stage of AAA and then modified HNA
opts.reg = 5e-10;
es = zeros(size(deg)); ts = zeros(3,length(deg));
for i = 1:length(deg)
    opts.tol = eps*1e4;
    opts.maxit = deg(i);
    [R,~,~,Ap,Bp,Cp,t1,H,Gstable,maxerr] = block_aaa(FF,samples,opts);
    Ga = @(z) Cp * ((((z-1)/(z+1)) * eye(length(Ap)) - Ap) \ Bp);
    Ap = blkdiag(A1,Ap); Bp = [B1;Bp]; Cp = [C1,Cp];
    [Ap,Bp,Cp,Dp,t2,t3] = modified_HNA(Ap,Bp,Cp,zeros(p,m),target_d,0.00001,1e-10);
    [Ap,Bp,Cp,Dp] = plane2disk(Ap,Bp,Cp,Dp);
    Gp = @(z) Dp + Cp * ((z*eye(length(Ap)) - Ap) \ Bp);
    es(i) = hankel_err(G,@(z) H((z-1)/(z+1))+Hinter((z-1)/(z+1))+Gp(z),err_conf);
    ts(1,i) = t1; ts(2,i) = t2; ts(3,i) = t3;
    fprintf(num2str(deg(i)));
end

%% Plot the results
if demo
    figure; hold on
    plot(deg,es,'-x','linewidth',2,'markersize',8,'displayname','error');
    plot([deg(1),deg(length(deg))],[err_bench,err_bench],'--','linewidth',2,...
        'displayname','theo. low. bd.');
    plot([deg(1),deg(length(deg))],[errb,errb],'--','linewidth',2,...
        'markersize',8,'displayname','BT error');
    xlabel('degree of AAA')
    ylabel('Hankel error')
    legend show
    grid on
    hold off

    figure; hold on
    plot(deg,ts(1,:),'-x','linewidth',2,'markersize',8,'displayname','t1');
    plot(deg,ts(2,:),'-x','linewidth',2,'markersize',8,'displayname','t2');
    plot(deg,ts(3,:),'-x','linewidth',2,'markersize',8,'displayname','t3');
    plot(deg,ts(1,:)+ts(2,:)+ts(3,:),'-x','linewidth',2,'markersize',8,...
        'displayname','total time');
    plot([deg(1),deg(length(deg))],[t_bench,t_bench],'--','linewidth',2,...
        'displayname','time of BT');
    xlabel('degree of AAA')
    ylabel('time')
    legend show
    grid on
    hold off
else
    fprintf('AAA + Glover error(s):');
    disp(es);
    fprintf('Theoretical lower bound:');
    disp(err_bench);
    fprintf('Balanced truncation error bound:');
    disp(errb);
    
    fprintf('\nAAA + Glover time:');
    disp(ts);
    fprintf('Balanced truncation time:');
    disp(t_bench);
end