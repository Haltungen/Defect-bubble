function A = MakeA_Bubble(omega,v0,vb,alp,delta,L,R,NN,N)
k0 = omega*v0;
kb = omega*vb;

M=2*NN+1;
A11c=zeros(M,M);
A12=zeros(M,M);
A21c=zeros(M,M);
A22=zeros(M,M);

A11d=zeros(M,M);
%A12d=zeros(M,M);
A21d=zeros(M,M);
%A22d=zeros(M,M);


%mm=1;


%%% store lattice sum data %%%



data_Sn_posi = zeros(2*NN+1,1);
data_Sn_nega = zeros(2*NN,1);
for j=1:2*NN+1
data_Sn_posi(j)=lattice_Sn(j-1,k0,alp,L,N);
end
for j=1:2*NN   
data_Sn_nega(j)=lattice_Sn(j,k0,[-alp(1), alp(2)],L,N);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=-NN:NN
    
    
    J0 = besselj(n, k0*R);
    H0 = besselh(n, 1, k0*R);
    
    Jb = besselj(n, kb*R);
    Hb = besselh(n, 1, kb*R);
    
    dH0 = 1/2*(besselh(n-1,1, k0*R) - besselh(n+1,1, k0*R));
    dJb = 1/2*(besselj(n-1, kb*R) - besselj(n+1, kb*R));
    
    
    c=1;
    
    %inside
    A12(NN+1+n,NN+1+n)=c*Jb*Hb;
    A22(NN+1+n,NN+1+n)=c*kb*dJb*Hb;
    
    %outside
    A11c(NN+1+n,NN+1+n)=c*J0*H0;
    A21c(NN+1+n,NN+1+n)=c*k0*J0*dH0;
    
    
    for m=-NN:NN
        
        %Snm0=makeSAn_new(n-m,k0,alp,L,N,mm);
        if (n-m) >= 0
        
        Snm0 = data_Sn_posi(n-m+1);
        else
            jj=-(n-m);
            Snm0 = data_Sn_nega(jj);
        end
        
        dJ0m = 1/2*(besselj(m-1,k0*R) - besselj(m+1,k0*R));
        J0m = besselj(m,k0*R);
        
        A11d(NN+1+m,NN+1+n)=c*J0*(-1)^(n-m)*Snm0*J0m;
        A21d(NN+1+m,NN+1+n)=c*k0*J0*(-1)^(n-m)*Snm0*dJ0m;
        
        %A21d(NN+1+m,NN+1+n)=c*k0*Snm0*dJ0m;
    end
        
       
        
end


%Ac=[A11c,A12c;A21c,A22c];
%Ad=[A11d,A12d;A21d,A22d];
%A=A21c+A21d;
A=[A11c+A11d,-A12;delta*(A21c+A21d),-A22];
end



