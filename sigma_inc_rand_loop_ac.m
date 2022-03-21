%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear;
format long

x = [1;0;0];
y = [0;1;0];
z = [0;0;1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sgrain=load('n100000-id1.ori');

K=1E5; % Number of materials
ValSigUni = 300; % Value of the uniaxial stress
angle = (-90:1:90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macroscopic stress tensor in the global frame (from FE in elastoplasticity)
SIGMA=zeros(3);
SIGMA(1,1)=ValSigUni;
% SIGMA(2,2)=ValSigUni;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic stiffnesses tensor 
% Ti beta--MPa
% Ravi from Bayesian inference 
C11=92600;
C12=82500;
C44=43500;
% Cu
% C11=170000;
% C12=124000;
% C44=75000;
% Al
% C11=107000; 
% C12=61000; 
% C44=28000;  

DLOCAL=zeros(6);
for j=1:3
    DLOCAL(j,j)=C11;
    for i=1:3
        if (i~=j)
            DLOCAL(i,j)=C12;
        end
    end
    DLOCAL(j+3,j+3)=C44;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GB volume fraction
% f=rand(K,1);
% suma=sum(f);
% for i=1:K
%     f(i)=f(i)/suma;
% end
f=zeros(K,1);
for i=1:K
    f(i)=1.0/K;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crossing matrix from crystals frame to the grain boundary frame
R=zeros(K,3,3);
T=zeros(K,3,3);

Nori=length(Sgrain);
for i=1:K
%     nb=ceil(rand*Nori);
    nb=i;
    phi1=Sgrain(nb,1)*pi/180;
    phi=Sgrain(nb,2)*pi/180;
    phi2=Sgrain(nb,3)*pi/180;
    R(i,:,:) = mrot(z,phi2)*mrot(x,phi)*mrot(z,phi1); % crossing matrix from the grain boundary frame to crystals frame
    T(i,:,:) = squeeze(R(i,:,:))^(-1); % crossing matrix from crystals frame to the grain boundary frame 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passing the 6x6 elastic stiffnesses matrix of crystals I, II, III and IV into the grain boundary frame
% conversion T1 (3x3) --> T1p (6x6)
Tp=zeros(K,6,6);

for k=1:K
    for j=1:3
        j1=1+floor(j/3);
        j2=2+floor(j/2);
        for i=1:3
            i1=1+floor(i/3);
            i2=2+floor(i/2);
            % convention 12 --> 4, 31 --> 5, 23 --> 6
            %         Tp(k,i,j)=T(k,i,j)^2;
            %         Tp(k,i,j+3)=2*T(k,i,j1)*T(k,i,j2);
            %         Tp(k,i+3,j)=T(k,i1,j)*T(k,i2,j);
            %         Tp(k,i+3,j+3)=T(k,i1,j1)*T(k,i2,j2)+T(k,i1,j2)*T(k,i2,j1);
            % convention 23 --> 4, 31 --> 5, 12 --> 6
            Tp(k,i,j)=T(k,i,j)^2;
            Tp(k,i,7-j)=2*T(k,i,j1)*T(k,i,j2);
            Tp(k,7-i,j)=T(k,i1,j)*T(k,i2,j);
            Tp(k,7-i,7-j)=T(k,i1,j1)*T(k,i2,j2)+T(k,i1,j2)*T(k,i2,j1);
        end
    end
end

C=zeros(K,6,6);
for m=1:K
    for i=1:6
        for j=1:6
            for k=1:6
                for l=1:6
                    C(m,i,j)=C(m,i,j)+DLOCAL(k,l)*Tp(m,i,k)*Tp(m,j,l);
                end
            end

        end
    end
end

% 4--6x6 elastic compliances matrix of crystals 
S=zeros(K,6,6);
for m=1:K
    S(m,:,:)=(squeeze(C(m,:,:)))^(-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncouche=ceil(log(K)/log(2))+1;
Kcouche=zeros(ncouche,1);
fcouche=zeros(ncouche,K);
Scouche=zeros(ncouche,K,6,6);
epcouche=zeros(ncouche,K,6);
epsilon_p=zeros(K,6);

KK=K;
ff=f;
Seff=S;
ep_eff=epsilon_p;

n=1;
Kcouche(n)=K;
fcouche(n,:,:)=f;
Scouche(n,:,:,:)=S;
epcouche(n,:,:)=epsilon_p;     

q=1;
if (K>1)
    while(q>0)
        [KK,ff,Seff,ep_eff,q] = Effprop(KK,ff,Seff,ep_eff);
        n=n+1;
        for i=1:KK
            Kcouche(n)=KK;
            fcouche(n,i)=ff(i);
            Scouche(n,i,:,:)=Seff(i,:,:);
            epcouche(n,i,:)=ep_eff(i,:);       
        end
    end
end

if (n ~= ncouche)
    disp('problem with ncouche')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VM=zeros(length(angle),1);
for ind=1:length(angle)
% Passing the macroscopic stress tensor into the grain boundary frame
    theta=angle(ind)*pi/180;
    Rot = mrot(z,theta);
    Sig = Rot*SIGMA*Rot^(-1);
 
% Macroscopic stress in vector notation
    Sig_Vect = [Sig(1,1);Sig(2,2);Sig(3,3);Sig(2,3);Sig(3,1);Sig(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stresses
    SIG=zeros(K,6);
    S2=squeeze(Scouche(ncouche,1,:,:));
    VonMises=zeros(K,1);
    SIG_SC=zeros(K,6);
    VonMises_SC=zeros(K,1);
    for i=1:K
        S1 = squeeze(S(i,:,:));
        GI=Gij(S1,S1,0);
        EP1 = squeeze(epsilon_p(i,:))';
        SIG(i,:) = Sig_Vect + GI*(S2-S1)*Sig_Vect;
        VonMises(i)=sqrt(0.5*( (SIG(i,1)-SIG(i,2))^2 + (SIG(i,2)-SIG(i,3))^2 + (SIG(i,3)-SIG(i,1))^2 ...
            + 6*(SIG(i,4)^2+SIG(i,5)^2+SIG(i,6)^2) ));
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VM(ind)=max(VonMises); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voigt_ind = [ 1 6 5
              6 2 4
              5 4 3 ];
          
C4=zeros(K,3,3,3,3);
for m=1:K
    for i=1:3, for j=1:3, for k=1:3, for l=1:3
           I = voigt_ind(i,j);
           J = voigt_ind(k,l);
           C4(m,i,j,k,l) = C(m,I,J);
    end; end; end; end
end

C_Voigt4 = zeros(3,3,3,3);
for m=1:K
    C_Voigt4 = C_Voigt4 + squeeze(C4(m,:,:,:,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ac = [1;2;5;10;50;100];

VM_SC=zeros(length(angle),length(ac));

for q=1:length(ac)
    C_eff = zeros(6,6);
    C_eff4 = C_Voigt4;
    for i=1:3, for j=1:3, for k=1:3, for l=1:3,
       I = voigt_ind(i,j);
       J = voigt_ind(k,l);
       C_eff(I,J) = C_eff4(i,j,k,l);
    end; end; end; end
    S_eff = inv(C_eff);

    B = zeros(K,6,6);
    testc = 1;
    while (testc > 1E-10)  
        [SE, P] = eshelby_tensor_aniso (C_eff4, ac(q), 1, ac(q), 1e-11);

        Mtild = (eye(6,6)-SE)\SE*S_eff;

        S_eff_old = S_eff;
        S_eff = zeros(6,6);
        acc = zeros(6,6);
        parfor m=1:K
            B(m,:,:) = (squeeze(S(m,:,:))+Mtild)\(S_eff_old+Mtild);
            S_eff = S_eff + f(m)*squeeze(S(m,:,:))*squeeze(B(m,:,:));
            acc = acc + f(m)*squeeze(B(m,:,:));
        end
        S_eff = S_eff / acc;

        testc = abs(norm(S_eff) - norm(S_eff_old))/norm(S_eff);

        C_eff = inv(S_eff);
        for i=1:3, for j=1:3, for k=1:3, for l=1:3,
           I = voigt_ind(i,j);
           J = voigt_ind(k,l);
           C_eff4(i,j,k,l)=C_eff(I,J);
        end; end; end; end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ind=1:length(angle)
% Passing the macroscopic stress tensor into the grain boundary frame
        theta=angle(ind)*pi/180;
        Rot = mrot(z,theta);
        Sig = Rot*SIGMA*Rot^(-1);

% Macroscopic stress in vector notation
        Sig_Vect = [Sig(1,1);Sig(2,2);Sig(3,3);Sig(2,3);Sig(3,1);Sig(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stresses
        SIG=zeros(K,6);
        VonMises_SC=zeros(K,1);
        for i=1:K
            SIG_SC(i,:) = squeeze(B(i,:,:))*Sig_Vect;
            VonMises_SC(i)=sqrt(0.5*( (SIG_SC(i,1)-SIG_SC(i,2))^2 + (SIG_SC(i,2)-SIG_SC(i,3))^2 + (SIG_SC(i,3)-SIG_SC(i,1))^2 ...
                + 6*(SIG_SC(i,4)^2+SIG_SC(i,5)^2+SIG_SC(i,6)^2) ));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VM_SC(ind,q)=max(VonMises_SC);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gca,'FontSize',20)
hold on
plot(angle,VM,'or','markersize',10,'Linewidth',1)
for q=1:length(ac)
    plot(angle,VM_SC(:,q),'s','markersize',10,'Linewidth',1)
end
xlabel('\theta (°)')
ylabel('Maximum von Mises stress (MPa)')
legend('off')
legend('Laminate','SC a/c=1','SC a/c=2','SC a/c=5','SC a/c=10','SC a/c=50','SC a/c=100')
legend('boxoff')
box on






