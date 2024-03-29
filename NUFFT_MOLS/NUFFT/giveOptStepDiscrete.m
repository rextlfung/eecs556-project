function [Kernel,calcweight,error,fn,step] = giveOptStepDiscrete(DFTMtx, fn,K,N,Ofactor, Olderror, Oldfn, a,Order)

for step = [0.5.^([0:1:30]),0],
    fn = step*a + (1-step)*Oldfn;
    [Kernel,calcweight,error] = calcKernelDiscretemod(DFTMtx,fn,K,N,Ofactor,Order);
    if(error<Olderror)
        break;
    end
end

