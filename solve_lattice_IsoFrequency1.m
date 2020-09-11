% solve_beam
clear;
tic;

%%

sample = 0.1;
u11 = (0.01:sample:pi)';
u22 = (0.01:sample:pi)';

index = 0;
u1 = [];
u2 = [];
for index1 = 1:length(u11)
    for index2 = 1:length(u22)
      index = index+1;
      u1(index) = u11(index1);
      u2(index) = u22(index2);
        
    end
end

feq = cell(length(u1),1);

%% Function Handler for spectral element
parfor iSample = 1:length(u1)

    funKd = @Main_function;
    fn   = 60;        % No. of samples in frequency domain
    fmin = 0.001;   % Minimum frequency (Hz)
    fmax = 1000;    % Maximum frequency (Hz)

    % interpolation points on circle
    radius = (fmax - fmin)/2;
    center = fmin + radius;
    z = exp(2i*pi*(0:fn-1)/(fn));
    f = radius*z + center;

    % samples of Kd
    Kd = funKd(f,u1(iSample),u2(iSample));
    s = size(Kd{1},1);


    % pole handling
    e = ones(s,1)/sqrt(s);
    ff = cellfun(@(x) e'*x*e,funKd(f,u1(iSample),u2(iSample)));
    [~,a,b,~,~,~,~] = ratdisk(ff,floor((fn-1)/2),floor((fn-1)/2),(fn-1),1e-14);
    poles = roots(b(end:-1:1));

    % remove poles
    P = Kd;
    for i = 1:length(P)
        for j = 1:length(poles)
            P{i} = P{i}*(z(i) - poles(j));
        end
        P{i} = P{i};
    end

    % scaling
    maxnrmP = max(cellfun(@(x) norm(x),P));
    for i = 1:length(P)
        P{i} = P{i}/maxnrmP;
    end


     % solve nlep
    [C0,C1] = linlagr(P,z,z);   % linearization
    [V,D] = eig(C0,C1);         % solve glep

    lam = radius*diag(D) + center;
    err = zeros(size(lam));

    for i = 1:length(lam)
        u = V(1:s,i);
        u = u/norm(u);
        Kdlam = funKd(lam(i),u1(iSample),u2(iSample));
        err(i) = norm(Kdlam{1}*u);
    end


    % output
    tol = 1e-6;
    sol{iSample} = lam(err<tol);


    % save solve_frame_linearization;
    % plot
%     figure
%     hold off; plot(radius*exp(2i*pi*(0:fn-1)/(fn)) + center,'-k');
%     hold on;
%     plot(lam,'*b');
%     plot(lam(err<tol),'*b');
%     plot(radius*poles + center+1i*1E-10,'or');
%     xlim([fmin fmax])
%     ylim([-fmax/2 fmax/2])

    %sol =   sol(unique(round(sort(real(sol)),4))>=0);
    feq{iSample} = unique(round(sort(real(sol{iSample})),4));
    feq{iSample} = feq{iSample}(feq{iSample}>=0);
    
    
end


save result_Honeycomb_ISO_0to1000.mat


for index = 1:size(feq,1)
    size_feq(index) = size(feq{index},1);
end
min(size_feq)

