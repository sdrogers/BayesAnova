function out = beta_clusterer(Q)


K = 10;
alpha = 1;
N = size(Q,1);
Z = rand(N,K);
Z = (Z==repmat(max(Z,[],2),1,K));


a_in = 10;b_in = 1;
a_out = 1;b_out = 10;

p0in = 0.01;
p0out = 0.5;


% Compute the pairwise likelihoods
Lin = lbetapdf(Q,a_in,b_in);
Lout = lbetapdf(Q,a_out,b_out);
Lin(Q==0) = log(p0in);
Lout(Q==0) = log(p0out);
NITS = 100;

baseLike = sum(Lout,2);
baseLike = baseLike - diag(Lout);

sZ = sum(Z,1);


for s = 1:NITS
    order = randperm(N);
    for n = 1:N
        this = order(n);
        thisk = find(Z(this,:));
        Z(this,thisk) = 0;
        sZ(thisk) = sZ(thisk) - 1;
        if sZ(thisk) == 0
            K = K - 1;
            Z(:,thisk) = [];
            sZ(thisk) = [];
        end
        
        prior = [sZ alpha];
        
        Like = repmat(baseLike(this),1,K);
        Like = Like - sum(Z.*repmat(Lout(:,this),1,K),1);
        Like = Like + sum(Z.*repmat(Lin(:,this),1,K),1);
        Post = [Like baseLike(this)] + log(prior);
        Post = exp(Post - max(Post));
        newk = find(rand<=cumsum(Post),1);
        
        
        if newk <= K
            Z(this,newk) = 1;
            sZ(newk) = sZ(newk) + 1;
        else
            Z(this,end+1) = 1;
            sZ(1,end+1) = 1;
            K = K + 1;
        end
        
        
        
        
    end
    
    [szI,p] = sort(sZ,'descend');
    I = [];
    for k = 1:K
        I = [I;find(Z(:,p(k)))];
    end
    imagesc(Q(I,I));drawnow
end



function l = lbetapdf(X,al,be)

l = gammaln(al+be) - gammaln(al) - gammaln(be) + ...
    (al-1).*log(X) + (be-1).*log(1-X);
