%% Do an individual ANOVA on each row in the file
clear all;
close all;

%% Load the data
dataname = '/Users/simonrogers/Dropbox/Meta_clustering/Human_PE/Positive_processed.csv';
X = importdata(dataname);

a = X.data(1,:); % Cases v controls
b = X.data(2,:); % 16 v 28 weeks

ma = max(X.data);
Y = X.data(3:end,:);
for i = 1:size(Y,2)
    temp = ~isnan(Y(:,i));
    s(i) = sum(Y(temp,i));
end
%%
Y = Y./repmat(s,size(Y,1),1);

% Y = X.data(3:end,:);
% ma = max(Y,[],1);
% Y = Y./repmat(ma,size(Y,1),1);

Y = log(Y);


pos = find(a==0 & b==0);
m = mean(Y(:,pos),2);
Y = Y - repmat(m,1,size(Y,2));
%%
NPeaks = size(Y,1);

all_alphasig = [];

for n = 1:NPeaks
    if ~any(isnan(Y(n,:)))
        out = oneDBayesAnova(Y(n,:),a,b,'plots',0);
        s = sprintf('Peak ID: %s, Alphasig: [%g,%g], ',X.textdata{n},out.alphasig)
        all_alphasig(n,:) = out.alphasig;
        all_betasig(n,:) = out.betasig;
        all_alphabetasig(n,:) = out.alphabetasig;
%         if any(out.alphasig<0.1)
%             keyboard
%         end
    else
        s = sprintf('Peak ID: %s, has nans',X.textdata{n});
        all_alphasig(n,:) = [nan nan];
        all_betasig(n,:) = [nan nan];
        all_alphabetasig(n,:) = [nan nan];
    end
end
    

%% Make a plot
pos = find(all_alphabetasig(:,1)<0.05);
for p = pos';
figure(2)
hold off
plot(Y(p,:),'ko');
hold all
out = oneDBayesAnova(Y(p,:),a,b,'plots',1,'gamma',10);
ma = mean(out.alpha_all);
mb = mean(out.beta_all);
mab = mean(out.alphabeta_all);
figure(2)
plot([1:72],a.*ma + b.*mb + a.*b.*mab,'r');
for i = 18.5:18:70
	plot([i i],ylim,'k--');
end
text(2,0,'16 week cases');
text(20,0,'16 week controls');
text(38,0,'28 week cases');
text(54,0,'28 week controls');
title(sprintf('Peak %s',X.textdata{p + 2}));
pause
end