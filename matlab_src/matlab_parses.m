%% Parse the data
clear all;
close all;
datapath = '/Users/simonrogers/Dropbox/Meta_clustering/Human_PE/';
fileList = importdata([datapath 'fileList.txt']);

nFiles = length(fileList);

X = [];

for n = 1:nFiles
    temp = importdata([datapath fileList{n}]);
    temp2 = full(sparse(temp.data(:,1),repmat(1,size(temp.data(:,1)),1),temp.data(:,4)));
    if length(temp2)<size(X,1)
        temp2 = [temp2;zeros(size(X,1)-length(temp2),1)];
    end
    X = [X temp2];
    a(n) = 0;
    q = strfind(fileList{n},'cases') ;
    if ~isempty(q)
        a(n) = 1
    end
    b(n) = 0;
    q = strfind(fileList{n},'28') ;
    if ~isempty(q)
        b(n) = 1;
    end
end

