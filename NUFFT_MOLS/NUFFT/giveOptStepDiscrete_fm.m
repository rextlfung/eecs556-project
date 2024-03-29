function [Kernel,calcweight,error,fn,step] = giveOptStepDiscrete_fm(DFTMtx, fn,K,N,Ofactor, Olderror, Oldfn, a,Order,H)

for step = [0.5.^([0:1:30]),0],
    fn = step*a + (1-step)*Oldfn;
    [Kernel,calcweight,error] = calcKernelDiscretemod_fm(DFTMtx,fn,K,N,Ofactor,Order,H);
    if(Olderror-error>0)
        break;
    end
end

