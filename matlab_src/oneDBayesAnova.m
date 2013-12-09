function out = oneDBayesAnova(Y,a,b,varargin)

% Normalise the data
pos = find(a==0 & b==0);
m = mean(Y(:,pos),2);
Y = Y - m;

gamma = 1;
NSAMPS = 1000;
NBURN = 500;
plots = 1;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'gamma'
            gamma = varargin{i+1};
        case 'nsamps'
            NSAMPS = varargin{i+1};
        case 'nburn'
            NBURN = varargin{i+1};
        case 'plots'
            plots = varargin{i+1};
    end
end


alpha = randn;
beta = randn;
alphabeta = randn;

out.alpha_all = zeros(NSAMPS-NBURN,1);
out.beta_all = zeros(NSAMPS-NBURN,1);
out.alphabeta_all = zeros(NSAMPS-NBURN,1);


for s = 1:NSAMPS

    qma = Y - (beta*b + alphabeta*(a.*b));
    ssa = 1./(1 + gamma*sum(a));
    mua = ssa*gamma*sum(qma.*a);

    alpha = randn.*sqrt(ssa) + mua;
    
    qmb = Y - (alpha*a + alphabeta*(a.*b));
    ssb = 1./(1 + gamma*sum(b));
    mub = ssb*gamma*sum(qmb.*b);

    beta = randn.*sqrt(ssb) + mub;
    
    
    qmab = Y - (alpha*a + beta*b);
    ssab = 1./(1 + gamma*sum(a.*b));
    muab = ssab*gamma*sum(qmab.*a.*b);

    alphabeta = randn.*sqrt(ssab) + muab;
    
    if(s>NBURN)
        out.alpha_all(s-NBURN) = alpha;
        out.beta_all(s-NBURN) = beta;
        out.alphabeta_all(s-NBURN) = alphabeta;
    end
end

if plots
    figure(1);hold off
    [f,xi] = ksdensity(out.alpha_all);
    plot(xi,f);
    hold all
    [f,xi] = ksdensity(out.beta_all);
    plot(xi,f);
    [f,xi] = ksdensity(out.alphabeta_all);
    plot(xi,f);
    xl = xlim;
    xr = xl(1):0.01:xl(2);
    plot(xr,normpdf(xr,0,1),'k--');
    legend('Alpha','Beta','AlphaBeta','Prior');
end

out.alphasig = [mean(out.alpha_all>0) mean(out.alpha_all<0)];
out.betasig = [mean(out.beta_all>0) mean(out.beta_all<0)];
out.alphabetasig = [mean(out.alphabeta_all>0) mean(out.alphabeta_all<0)];