function [velL,velR]=profileGenAsym(mainSpeeds,perturbationSizes,pertDuration,washDuration, reps)
%v1: 21/9/2015 - Pablo
%v2: 16/8/2016 - Pablo
%% Initialize blocks:
M=length(mainSpeeds);
N=length(perturbationSizes);
P=pertDuration;
W=washDuration;
R=reps;
velL=nan(P+W,2*M*N,R); %trial duration x number of trials x number of blocks
velR=nan(P+W,2*M*N,R);

%% create the mother block:
blockL=nan(P+W,M,N,2);
blockR=nan(P+W,M,N,2);

for i=1:M
    for j=1:N
        for k=1:2 %L and R
            columnP=mainSpeeds(i)*ones(P+W,1);
            columnC=mainSpeeds(i)*ones(P+W,1) + [zeros(W,1); perturbationSizes(j)*ones(P,1)];
            switch k
                case 1
                    blockL(:,i,j,k)=columnP;
                    blockR(:,i,j,k)=columnC;
                case 2
                    blockR(:,i,j,k)=columnP;
                    blockL(:,i,j,k)=columnC;
            end
        end
    end
end

%% Permute order in each repetition:
for i=1:R %For each repetition
    KK=randperm(M*N*2);
    velL(:,:,i)=blockL(:,KK);
    velR(:,:,i)=blockR(:,KK);
end

%% 
%Set self-controlled strides to NaN, as required by controller
velL(W+4:end,:,:)=NaN;
velR(W+4:end,:,:)=NaN;

%%
velL=[velL(:); mainSpeeds(1)*ones(10,1)];
velR=[velR(:); mainSpeeds(1)*ones(10,1)];
end