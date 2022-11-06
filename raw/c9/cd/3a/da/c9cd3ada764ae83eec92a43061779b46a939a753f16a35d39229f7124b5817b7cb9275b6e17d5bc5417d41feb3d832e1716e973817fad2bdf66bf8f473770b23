%RWG33 FREQUENCY LOOP FOR THE SCATTERING PROBLEM 
%   Calculates the voltage in the feed of the receiving 
%   antenna if the incident signal (E-field) is given.
%   Performs calculations at all required frequencies
%
%   Uses mesh2.mat (created by RWG2), current.mat (created by RWG31), 
%   and radiatedfield.mat (created by RWG32) as inputs
% 
%   The following parameters need to be specified prior to 
%   calculations:
%   
%   Number of frequency steps       NumberOfSteps
%   Lower frequency                 FreqStart
%   Upper frequency                 FreqStop
%   Dielectric constant (SI)        epsilon_
%   Magnetic permeability (SI)      mu_
%   Wave vector                     kv
%
%   Copyright 2002 AEMM. Revision 2002/03/26 
%   Chapter 8

%load the data
load('mesh2');
load('current')
load('radiatedfield');

%Frequency series parameters
NumberOfSteps=500;
FreqStart   =12.5e6;  %in Hz
FreqStop    =6250e6;  %in Hz
step=(FreqStop-FreqStart)/(NumberOfSteps-1);

%EM parameters
epsilon_    =8.854e-012;
mu_         =1.257e-006;
%Speed of light 
c_=1/sqrt(epsilon_*mu_);
%Free-space impedance 
eta_=sqrt(mu_/epsilon_);

%Contemporary variables - metal impedance matrix
for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   %[3 9 EdgesTotal]
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  %[3 9 EdgesTotal]
end

%Frequency series
T0=cputime;
for FF=1:NumberOfSteps    
    FF
    f(FF)       =FreqStart+step*(FF-1);
    omega       =2*pi*f(FF);
    k           =omega/c_;
    K           =j*k;
  
    Constant1   =mu_/(4*pi);
    Constant2   =1/(j*4*pi*omega*epsilon_);
    Factor      =1/9;    
    FactorA     =Factor*(j*omega*EdgeLength/4)*Constant1;
    FactorFi    =Factor*EdgeLength*Constant2;
    FactorA     =FactorA.';
    FactorFi    =FactorFi.';
    
    Z   =  impmet( EdgesTotal,TrianglesTotal,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);   
    
   
    d=[0 0 -1];     
    kv=k*d;
    Pol=E(1:3,FF).';
    
    for m=1:EdgesTotal    
        ScalarProduct   =sum(kv.*Center(:,TrianglePlus(m))');
        EmPlus          =Pol.'*exp(-j*ScalarProduct);      
        ScalarProduct   =sum(kv.*Center(:,TriangleMinus(m))');
        EmMinus         =Pol.'*exp(-j*ScalarProduct);      
        ScalarPlus      =sum(EmPlus.* RHO_Plus(:,m));
        ScalarMinus     =sum(EmMinus.*RHO_Minus(:,m));
        V(m)            =EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2);   
    end
        
    %Solution of MoM equations
    I=Z\V.';
    
    %Find the antenna feed position (receiving antenna)
    FeedPoint=[0; 0; 0];
    for m=1:EdgesTotal
        V(m)=0;
        Distance(:,m)=0.5*sum(p(:,Edge_(:,m)),2)-FeedPoint;
    end
    [Y,INDEX]=sort(sum(Distance.*Distance));
    Index=INDEX(1);                 %Center feed - dipole
    %Index=INDEX(1:2);              %Probe feed - monopole
      
    %Find the voltage across the antenna feed of the receiving
    %antenna
    Imp=Impedance(FF);
    FeedCurReceived     =I(Index)*EdgeLength(Index);
    FeedVolReceived     =FeedCurReceived*Imp;
    OutputVoltage(FF)   =FeedVolReceived;
    PowerConjMatch(FF)  =1/8*(FeedVolReceived*conj(FeedVolReceived))/real(Imp);    
    
    T=cputime-T0
end

%Save result
FileName='receivedfield.mat'; 
save(FileName, 'f','NumberOfSteps','FreqStart','FreqStop','OutputVoltage','PowerConjMatch',...
    'Impedance','Index');
plot(f,PowerConjMatch);

            