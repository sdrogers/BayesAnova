%% Test the two-way model model
clear all;
close all;

%% Generate some data
K = 6; % Number of clusters
M = 40; % Number of samples
a = [repmat(0,M/2,1);repmat(1,M/2,1)];
b = [repmat(0,M/4,1);repmat(1,M/4,1)];
b = repmat(b,2,1);

low_gamma = 1;
high_gamma = 10;
cov_gamma = 0.1;

true_alpha = [5;0;0;-2;3;0];
true_beta = [0;5;3;4;0;0];
true_alphabeta = zeros(K,1);
true_alphabeta(end) = -3;

true_mu = zeros(K,1);
true_x = zeros(K,M);
% Generate the true x
for k = 1:K
    true_x(k,:) = randn(1,M)./sqrt(low_gamma) + true_mu(k) + true_alpha(k).*a' + true_beta(k).*b' + true_alphabeta(k).*a'.*b';
end


% Generate the data
Y = [];
flatZ = [];
for k = 1:K
    Nk(k) = poissrnd(10);
    Y = [Y;randn(Nk(k),M)./sqrt(high_gamma) + repmat(true_x(k,:),Nk(k),1)];
    flatZ = [flatZ;repmat(k,Nk(k),1)];
end
N = size(Y,1);
true_Z = full(sparse([1:N]',flatZ,1));

%% Run the sampler

alpha = randn(K,1);
beta = randn(K,1);
alphabeta = zeros(K,1); % Ignore this for now
X = zeros(K,M);

NSAMPS = 1000;

alpha_all = zeros(NSAMPS,K);
beta_all = zeros(NSAMPS,K);
alphabeta_all = zeros(NSAMPS,K);

Z = rand(N,K);
Z = (Z==repmat(max(Z,[],2),1,K));

for s = 1:NSAMPS
    % Update x
    for k = 1:K
        for m = 1:M
            q = a(m)*alpha(k) + b(m)*beta(k) + a(m)*b(m)*alphabeta(k);
           ssx = 1/(1 + high_gamma*sum(Z(:,k)));
           mux = ssx*(q + high_gamma*sum(Z(:,k).*Y(:,m)));
           X(k,m) = randn.*sqrt(ssx) + mux;
        end
    end
    
    % Update alpha
    qma = X - (beta*b' + alphabeta*(a.*b)');
    ssa = 1./(1 + low_gamma*sum(a));
    mua = ssa*low_gamma*sum(qma.*repmat(a',K,1),2);
    alpha = randn(K,1).*sqrt(ssa) + mua;
   
%     for k = 1:K
%         qma = X(k,:) - (b'.*beta(k) + a'.*b'.*alphabeta(k));
%         ssa = 1/(1 + low_gamma*sum(a));
%         mua = ssa*low_gamma*sum(a'.*qma);
%         alpha(k) = randn.*sqrt(ssa) + mua;
%     end
    
    qmb = X - (alpha*a' + alphabeta*(a.*b)');
    ssb = 1./(1 + low_gamma*sum(b));
    mub = ssb*low_gamma*sum(qmb.*repmat(b',K,1),2);
    beta = randn(K,1).*sqrt(ssb) + mub;


%     % Update beta
%     for k = 1:K
%         qmb = X(k,:) - (a'.*alpha(k) + a'.*b'.*alphabeta(k));
%         ssb = 1/(1 + low_gamma*sum(b));
%         mub = ssb*low_gamma*sum(b'.*qmb);
%         beta(k) = randn.*sqrt(ssb) + mub;        
%     end
    

    qmab = X - (beta*b' + alpha*a');
    ssab = 1./(1 + low_gamma*sum(a.*b));
    muab = ssab*low_gamma*sum(qmab.*repmat((a.*b)',K,1),2);
    alphabeta = randn(K,1).*sqrt(ssab) + muab;


%     % Update alphabeta
%     for k = 1:K
%         qmab = X(k,:) - (a'.*alpha(k) + b'.*beta(k));
%         ssab = 1/(1 + low_gamma*sum(a.*b));
%         muab = ssab*low_gamma*sum((a.*b)'.*qmab);
%         alphabeta(k) = randn.*sqrt(ssab) + muab;        
%     end
    
    
    alpha_all(s,:) = alpha';
    beta_all(s,:) = beta';
    alphabeta_all(s,:) = alphabeta';
    
    % Update the clustering matrix
    Like = zeros(N,K);
    for k = 1:K
        Like(:,k) = lnpdf(Y,X(k,:),high_gamma);
    end
    Probs = exp(Like - repmat(max(Like,[],2),1,K));
    for n = 1:N
        newk = find(rand<=cumsum(Probs(n,:)),1);
        Z(n,:) = 0;
        Z(n,newk) = 1;
    end
    
    imagesc(Z);drawnow
    
end

figure(1);hold off
legcell = {};
for k = 1:K
    [f,xi] = ksdensity(alpha_all(:,k));
    plot(xi,f);
    hold all
    legcell{k} = num2str(true_alpha(k));
end
legend(legcell)
title('Alpha');

figure(2);hold off
legcell = {};
for k = 1:K
    [f,xi] = ksdensity(beta_all(:,k));
    plot(xi,f);
    hold all
    legcell{k} = num2str(true_beta(k));
end
legend(legcell)
title('Beta');


figure(3);hold off
legcell = {};
for k = 1:K
    [f,xi] = ksdensity(alphabeta_all(:,k));
    plot(xi,f);
    hold all
    legcell{k} = num2str(true_alphabeta(k));
end
legend(legcell)
title('AlphaBeta');

