function [walkPaths,endStep,endSign] = simulateRandomWalk(threshold, driftRate, noiseRate, simTime, NstepsPerTimeUnit,Niter,bias,fixedNoise)

T=simTime;
dt=1/(NstepsPerTimeUnit); %Integration time-step
Nsteps=ceil(T/dt);
mu=driftRate;
sigma=noiseRate;

    endStep=nan(Niter,1);
    endSign=nan(Niter,1);
    walkPaths=nan(Niter,Nsteps+1);
    for iter=1:Niter
        w=nan(Nsteps+1,1);
        w(1)=bias;
        for k=[1:Nsteps] %Time-steps
            w(k+1)=w(k)+((sqrt(dt)*sigma)*randn(1)+mu*dt); 
            %Multyplying by dt so that mu and sigma are inforamtion/uncertainty per unit of time, as opposed to per step
            %In this way, the limit when Nstep is large is a point
            %process. Otherwise, the limit would be a trivial case in
            %which the walk ends at k=1 always and where the
            %probability of a '->' choice is the right-tail of a normal
            %distribution:
            %.5*(1-erf((th-mu)/(sqrt(2)*sigma))
            if abs(w(k+1)+fixedNoise*randn(1))>(threshold)
                endStep(iter)=k*dt; %Not allowing subjects to respond in the very first 2% of the time interval
                endSign(iter)=sign(w(k+1));
                break
            end
        end
        walkPaths(iter,:)=w;
    end

end

