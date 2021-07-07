function y1=getData2Mtrx(Data,segmNames1,variable1,variable2,stat1)

y1=[];
for n=1:length(segmNames1)
y1(n,:)=[Data.(segmNames1{n}).(variable1).(variable2).(stat1)];
end