% Perturbation size responses for different time bins

figure;
aux=unique((allAuxZ));
hold on; 
h=histogram(allAuxZ(allAuxX<3 & allAuxX>2),[-400;aux(2:2:end)+5;400],'EdgeColor','none','Normalization','probability');
histogram(allAuxZ(allAuxX<4 & allAuxX>3),h.BinEdges,'EdgeColor','none','Normalization','probability')
histogram(allAuxZ(allAuxX<6 & allAuxX>5),h.BinEdges,'EdgeColor','none','Normalization','probability')
plot(aux,0,'r*')

legend('2-3s','3-4s','5-6s')

%%

aux=unique(abs(allAuxZ));
figure;
hold on; 
h=histogram(abs(allAuxZ(allAuxX<3 & allAuxX>2)),[0;aux(1:1:end)+5;aux(end)+5],'EdgeColor','none','Normalization','probability');
histogram(abs(allAuxZ(allAuxX<4 & allAuxX>3)),h.BinEdges,'EdgeColor','none','Normalization','probability')
histogram(abs(allAuxZ(allAuxX<6 & allAuxX>5)),h.BinEdges,'EdgeColor','none','Normalization','probability')
plot(aux,0,'r*')
axis tight
legend('2-3s','3-4s','5-6s')