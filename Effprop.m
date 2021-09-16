function [KK,ff,Seff,ep_eff,q] = Effprop(K,f,S,ep)

Seff=zeros(K,6,6);
ep_eff=zeros(K,6);
ff=zeros(K,6);
KK=0;
q=1;
for i=1:K-1
    if (mod(i,2)==1)
        KK=KK+1;
        fv=f(i)/(f(i)+f(i+1));
        S1=squeeze(S(i,:,:));
        S2=squeeze(S(i+1,:,:));
        ep1=squeeze(ep(i,:))';
        ep2=squeeze(ep(i+1,:))';
        G=Gij(S1,S2,fv);
        Seff(KK,:,:) = fv*S1+(1-fv)*S2 - fv*(1-fv)*(S2-S1)*G*(S2-S1);
        ep_eff(KK,:) = fv*ep1+(1-fv)*ep2 - fv*(1-fv)*(S2-S1)*G*(ep2-ep1);
        ff(KK)=f(i)+f(i+1);
    end
end
if (mod(K-1,2)==0)
    KK=KK+1;
    Seff(KK,:,:)=S(K,:,:);
    ep_eff(KK,:)=ep(K,:);
    ff(KK)=f(K);
end
if(KK==1)
    q=0;
end
