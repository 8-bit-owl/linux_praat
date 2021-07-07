function dataout=findclusters(data)
%FINDCLUSTERS get 2 clusters and sort them
%
%
%
% ejh 2019


obj = fitgmdist(data,2);
idx = cluster(obj,data);

cluster1 = data(idx == 1);
cluster2 = data(idx == 2);
gcluster1 = fitdist(cluster1,'Normal');
gcluster2 = fitdist(cluster2,'Normal');
cluster1stats= basicstats(cluster1);
cluster2stats= basicstats(cluster2);
if cluster1stats.mean>cluster2stats.mean % find the smaller mean
    dataout.low.cluster= cluster2;
    dataout.low.gcluster= gcluster2;
    dataout.low.stats= cluster2stats;
    dataout.high.cluster= cluster1;
    dataout.high.gcluster= gcluster1;
    dataout.high.stats=cluster1stats;
else
    dataout.low.cluster= cluster1;
    dataout.low.gcluster= gcluster1;
    dataout.low.stats= cluster1stats;
    dataout.high.cluster= cluster2;
    dataout.high.gcluster= gcluster2;
    dataout.high.stats=cluster2stats;
end
