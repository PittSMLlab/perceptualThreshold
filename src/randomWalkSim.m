%% TEst random walk integrate and fire

%%
th=1;
T=1; %End time, in a.u.
Nsteps=250;
dt=T/(Nsteps+1); %Integration time-step
M=20;
mus=[-1:.01:1]*th*M; 
sigmas=[3e-1,1.5,3]*th*M; %std of random walk step, per unit of time
legSigmas={};
for i=1:length(sigmas)
    legSigmas{i}=['\sigma=' num2str(sigmas(i),2)];
end
Niter=1e2;
allEndTime=nan(length(mus),length(sigmas),Niter);
allEndSign=nan(length(mus),length(sigmas),Niter);


for i=1:length(mus)
    mu=mus(i);
    for j=1:length(sigmas)
        sigma=sigmas(j);
        endCounter=nan(Niter,1);
        endSign=nan(Niter,1);
        for iter=1:Niter
            w=0;
            for k=[1:Nsteps] %Time-steps
                w=w+(sigma*randn(1)+mu)*dt; 
                %Multyplying by dt so that mu and sigma are inforamtion/uncertainty per unit of time, as opposed to per step
                %In this way, the limit when Nstep is large is a point
                %process. Otherwise, the limit would be a trivial case in
                %which the walk ends at k=1 always and where the
                %probability of a '->' choice is the right-tail of a normal
                %distribution:
                %.5*(1-erf((th-mu)/(sqrt(2)*sigma))
                if abs(w)>(th)
                    endCounter(iter)=max([k*dt,.02]); %Not allowing subjects to respond in the very first 2% of the time interval
                    endSign(iter)=sign(w);
                    break
                end
            end
        end
        allEndTime(i,j,:)=deal(endCounter);
        allEndSign(i,j,:)=deal(endSign);
    end
end

%% Compute some stats:
meanTime=nanmean(allEndTime,3);
stdTime=nanstd(allEndTime,[],3);
proportionRight=nanmean(allEndSign,3)/2 +.5;
stdRight=nanstd(allEndSign,[],3)/2+.5;
proportionNR=mean(isnan(allEndSign),3);
stdNR=std(isnan(allEndSign),[],3);
%% Do some plots
figure
subplot(2,2,1)
hold on
plot(mus,meanTime,'LineWidth',2)
cc=get(gca,'ColorOrder');
for i=1:length(sigmas)
   patch([mus'; mus(end:-1:1)'],[meanTime(:,i); meanTime(end:-1:1,i)]+[stdTime(:,i); -stdTime(end:-1:1,i)],cc(i,:),'EdgeColor','none','FaceAlpha',.3) 
end
legend(legSigmas)
xlabel('Speed diff (a.u.)')
ylabel('Mean response time')
set(gca,'YScale','log')

subplot(2,2,2)
hold on
plot(mus,proportionRight)
for i=1:length(sigmas)
    %p=fitPsychoGauss(mus,proportionRight(:,i),'MLE');
    %plot(mus,psychoGauss(p,mus),'k')
    %text(-M*.8, .8 -.1*i,['\mu=' num2str(p(1),2)],'Color',cc(i,:))
    %text(-M*.8, .75-.1*i,['\sigma_A=' num2str(p(2),2) '=' num2str(p(2)/sigmas(i)^2,2) '\sigma^2'],'Color',cc(i,:))
    p=fitPsycho(mus,proportionRight(:,i),'MLE');
    plot(mus,psycho(p,mus),'k')
    text(M*.1, .6 -.1*i,['a=' num2str(p(1),2)],'Color',cc(i,:))
    text(M*.1, .55-.1*i,['b=' num2str(p(2),2) '=' num2str(p(2)/sigmas(i)^2,2) '\sigma^2'],'Color',cc(i,:))

end
text(M*.1,.1,['Model: p=(1+e^{(x-a)/b})^{-1}'])
ylabel('% of rightward choices')
xlabel('Speed diff (a.u.)')

subplot(2,2,3)
hold on
plot(mus,proportionNR)
ylabel('% of NR')
xlabel('Speed diff (a.u.)')