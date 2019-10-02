%
%Null hypothesis: Choice probability depends on belt-speed difference only,
%and not on belt-speed ratio.
%Alternative: Anything else. Specifically, we expect that Weber's law is true, 
%and choice probabilities scale with the ratio of belt-speeds.

%Model of choice:
%p=1./(1+exp(s))
%Under the null, s=k.Dv
%Under the alt, s=k.v1/v2

%We collect 10 subjects in each of two groups. The speeds of group two are
%scaled by a factor of f with respect to the tested speeds of group one.
%For each subject, we test 13 different values of Dv (with the same v1+v2),
%repeating each 7 times.

f=1.2;
DvA=[-30,-150,-100,-75,-50,-25,0,25,50,75,100,150,300];
vA=1050;
DvB=DvA*f;
vB=vA*f;
Nsub=10;
Nreps=7;

%Simulate group 1:
k=15;
pA=1./(1+exp(-k*DvA/vA)); %k is set such that responses for 75 mm/s are about 75%, which is what we've measured
rA=binornd(Nsub*Nreps,pA);

%Simulate group 2 under null (only Dv matters):
pBnull=1./(1+exp(-k*DvB/vA));
rBnull=binornd(Nsub*Nreps,pBnull);

%Simualte group 2 under alt (Weber's law):
pBalt=pA;
rBalt=binornd(Nsub*Nreps,pBalt);


%% Do regressions to find best-fitting logistic:
rr=[rA;rBalt];
dv=[DvA;DvB];
for j=1:2
    responses=zeros(Nsub*Nreps,length(pA));
    for k=1:length(pA)
    responses(1:rr(j,k),k)=1;
    end
    aux=repmat(dv(j,:),Nsub*Nreps,1);
    aux=aux(:);
    responses=responses(:);
    tbl=table(responses,aux);
    mdl{j}=fitglm(tbl,'responses~aux-1','Distribution','binomial');
end
mdl{1}
mdl{2}

%NOTE: all these results are assumming that subjects have 0 bias, and it
%does not need to be estimated.

%TBD: repeat the process multiple times, decide on proper test to compare
%the two regression coefficient estimates, and decide on minimum value of f
%such that we can find a difference at least 80% of the time when testing
%for p<.05. This will ensure the 80% power.
%Then, repeat the process under the null, and check that only 5% of the
%time would we declare a difference when none exists.

%TBD2: this analysis treats all subjects as identical, and performs a
%group-level analysis. To add realism, we should model each subject's k as
%drawn from some distribution (normal?) and repeat the process. This should
%make things more difficult, as the distribution across subjects will
%probably make distributions wider and hence harder to distinguish.
%Additionally, I do not know how to properly do the estimation: we would
%need to consider a logistic model where subject ID is a (random) factor
%that changes the scaling with Dv (i.e. interaction term of ID and Dv), but
%the coefficient we care about is the average interaction, or something
%like that. Matlab will return the first subject's interaction as the main
%term, and all other subjects as departures with respect to the first
%subject. Not sure how to test with that data.
