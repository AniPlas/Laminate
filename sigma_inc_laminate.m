% This MATLAB code computes the effective behavior (effective stiffness tensor and the effective plastic strain tensor) and 
% the stress fields in an N-phase laminate with planar interface. 
% The laminate is submitted to a macroscopic homogeneous and remotely applied stress as well as to piecewise uniform plastic strains. 
% The phases of the laminate can have different grain volume fractions and can correspond to different materials or crystals
% The phases are assumed perfectly bonded  along the planar interface whose normal is along e2. 
% The code below is an application for a beta-titanium alloy with elongated grains loading in elasticity only, 
% modeled as laminate made of 100,000 different orientations and equal volume fraction.
% At the end, the maximal von Mises stress is plotted with respect to theta, a rotation around e3 of the uniaxial stress of magnitude 300 MPa.
% A reference for the model can be found in:

% The contracted Voigt notation (11:1, 22: 2, 33: 3, 23: 4, 31:5, 12: 6) is used for stresses and strains in vector notation and 
% stiffnesses and compliances in 6x6 matrix notation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear;
format long

x = [1;0;0];
y = [0;1;0];
z = [0;0;1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sgrain=load('n100000-cubic-uniform.ori');

K=1E5; % Number of materials
ValSigUni = 300; % Value of the uniaxial stress
angle = (-90:1:90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macroscopic stress tensor in the global frame 
SIGMA=zeros(3);
SIGMA(1,1)=ValSigUni;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic stiffnesses tensor 
% Ti beta--MPa
C11=92600;
C12=82500;
C44=43500;

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
% GB volume fraction
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
    TT=squeeze(Tp(m,:,:));
    C(m,:,:)=TT*DLOCAL*TT';
end


% 4--6x6 elastic compliances matrix of crystals 
S=zeros(K,6,6);
for m=1:K
    S(m,:,:)=(squeeze(C(m,:,:)))^(-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No plastic strain considered in this example
epsilon_p=zeros(K,6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of the laminate effective behavior from a multiple steps procedure

nstep=ceil(log(K)/log(2))+1; % total number of steps
Kstep=zeros(nstep,1); % number of phases at each step
fstep=zeros(nstep,K); % phase volume fraction at each step
Sstep=zeros(nstep,K,6,6); % phase elastic compliance tensor at each step
epstep=zeros(nstep,K,6); % phase plastic strain tensor at each step

% Initialization of phase properties at step n=1
n=1;
Kstep(n)=K;
fstep(n,:,:)=f;
Sstep(n,:,:,:)=S;
epstep(n,:,:)=epsilon_p;
KK=K;
ff=f;
Seff=S;
ep_eff=epsilon_p;
    
% Iterative loop: the effective behavior is known at the final nstep and Kstep = 1
q=1;
if (K>1)
    while(q>0)
        [KK,ff,Seff,ep_eff,q] = Effprop(KK,ff,Seff,ep_eff); % computation of effective properties of 2-phase laminate
        n=n+1;
        for i=1:KK
            Kstep(n)=KK;
            fstep(n,i)=ff(i);
            Sstep(n,i,:,:)=Seff(i,:,:);
            epstep(n,i,:)=ep_eff(i,:);       
        end
    end
end

if (n ~= nstep)
    disp('problem with nstep')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VM=zeros(length(angle),1); % maximum von Mises stress in the laminate

for ind=1:length(angle)
% Passing the macroscopic stress tensor into the grain boundary frame
    theta=angle(ind)*pi/180; % rotation angle around e3
    Rot = mrot(z,theta); % rotation matrix
    Sig = Rot*SIGMA*Rot^(-1); % macroscopic stress tensor
 
% Macroscopic stress in vector notation
    Sig_Vect = [Sig(1,1);Sig(2,2);Sig(3,3);Sig(2,3);Sig(3,1);Sig(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stresses
    SIG=zeros(K,6);
    S2=squeeze(Sstep(nstep,1,:,:)); % effective compliance matrix
	EP2=squeeze(epstep(nstep,1,:,:)); % effective plastic strain vector
    VonMises=zeros(K,1);
    for i=1:K
        S1 = squeeze(S(i,:,:));
        GI=Gij(S1,S2,0);
        EP1 = squeeze(epsilon_p(i,:))';
		SIG(i,:) = Sig_Vect + GI*( (S2-S1)*Sig_Vect + EP2-EP1 );  % computation of phase stress vector
        VonMises(i)=sqrt(0.5*( (SIG(i,1)-SIG(i,2))^2 + (SIG(i,2)-SIG(i,3))^2 + (SIG(i,3)-SIG(i,1))^2 ...
            + 6*(SIG(i,4)^2+SIG(i,5)^2+SIG(i,6)^2) )); % phase von Mises stress
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VM(ind)=max(VonMises); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gca,'FontSize',20)
hold on
plot(angle,VM,'sb','markersize',10,'Linewidth',3)
xlabel('\theta (°)')
ylabel('Maximum von Mises stress (MPa)')
legend('off')
box on







