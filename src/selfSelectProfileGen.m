p=cell(6,3);
p(:,1)={ones(1,20)};
p(:,3)={nan(1,20)};
p{1,2}=1.025*ones(1,2);
p{2,2}=1.1*ones(1,2);
p{3,2}=1*ones(1,2);
p{4,2}=1.075*ones(1,2);
p{5,2}=1.05*ones(1,2);
p{6,2}=1*ones(1,2);
p=p';
p=cell2mat(p(:)');
velR=p;
velL=2-p;
velL=repmat(velL,1,10); %10 reps
velR=repmat(velR,1,10);
save test7.mat velL velR

%%
Nreps=6;
p=cell(10,3);
p(:,1)={ones(1,12)}; %12 tied-belt strides, which include 2 for smooth return to base
p(:,3)={nan(1,15)}; %15 self-controlled strides
p{1,2}=1.1*ones(1,2);
p{2,2}=.975*ones(1,2);
p{3,2}=.925*ones(1,2);
p{4,2}=1.05*ones(1,2);
p{5,2}=1*ones(1,2);
p{6,2}=.9*ones(1,2);
p{7,2}=1.025*ones(1,2);
p{8,2}=1.075*ones(1,2);
p{9,2}=.95*ones(1,2);
p{10,2}=1*ones(1,2);
p=p';
p=cell2mat(p(:)');
velR=p;
velL=2-p;
velL=repmat(velL,1,Nreps); %6 reps
velR=repmat(velR,1,Nreps);
velL=[velL(end-28:end) velL]; %Adding null trial at the very beginning
velR=[velR(end-28:end) velR];
save test8.mat velL velR

%%
Nreps=6;
p=cell(10,3);
p(:,1)={ones(1,12)}; %12 tied-belt strides, which include 2 for smooth return to base
p(:,3)={nan(1,15)}; %15 self-controlled strides
p{1,2}=1.1*ones(1,2);
p{2,2}=.975*ones(1,2);
p{3,2}=.925*ones(1,2);
p{4,2}=1.05*ones(1,2);
p{5,2}=1*ones(1,2);
p{6,2}=.9*ones(1,2);
p{7,2}=1.025*ones(1,2);
p{8,2}=1.075*ones(1,2);
p{9,2}=.95*ones(1,2);
p{10,2}=1*ones(1,2);
p=p';
p=cell2mat(p(:)');
velR=p;
velL=2-p;
velL=repmat(velL,1,Nreps); %6 reps
velR=repmat(velR,1,Nreps);
velL=[velL(end-28:end) velL ones(1,30)]; %Adding null trial at the very beginning, and final washout
velR=[velR(end-28:end) velR ones(1,30)];
save test9.mat velL velR

%% Pilot 
Nreps=5;
p=cell(12,3);
p(:,1)={ones(1,12)}; %12 tied-belt strides, which include 2 for smooth return to base
p(:,3)={nan(1,15)}; %15 self-controlled strides
p{1,2}=.985*ones(1,2);
p{2,2}=1.1*ones(1,2);
p{3,2}=.975*ones(1,2);
p{4,2}=.925*ones(1,2);
p{5,2}=1.05*ones(1,2);
p{6,2}=1*ones(1,2);
p{7,2}=1.015*ones(1,2);
p{8,2}=.9*ones(1,2);
p{9,2}=1.025*ones(1,2);
p{10,2}=1.075*ones(1,2);
p{11,2}=.95*ones(1,2);
p{12,2}=1*ones(1,2);
p=p';
p=cell2mat(p(:)');
velR=p;
velL=2-p;
velL=repmat(velL,1,Nreps); %6 reps
velR=repmat(velR,1,Nreps);
velL=[velL(end-28:end) velL ones(1,30)]; %Adding null trial at the very beginning, and final washout
velR=[velR(end-28:end) velR ones(1,30)];
save thPilot_main.mat velL velR

clear velL velR
velL=ones(1,300);
velR=velL;
save thPilot_base.mat velL velR

clear velL velR
velL=[ones(1,30) nan(1,50) ones(1,12) nan(1,30) ones(1,12) nan(1,15) ones(1,12)];
velR=velL;
save thPilot_fam.mat velL velR

