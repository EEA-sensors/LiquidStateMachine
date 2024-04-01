%% Generate 2D spatially adjacent connected network (quad type)
% Net Size 2 element vector, W0 = 2x2 matrix Wt[EE EI;IE II]
function [X,Xn,W,R,E] = get3Dgrid(resSize,W0,fi,UC)
N = prod(resSize);
[RX,RY,RZ]= ndgrid(1:resSize(1),1:resSize(2),1:resSize(3));
%R(:,:,1) = RX; R(:,:,2) = RY; R(:,:,3) = RZ;
R=[]; R(:,1)=RX(:);R(:,2)=RY(:);R(:,3)=RZ(:);
R = reshape(R,numel(R)/3,3);
if(nargin<4 || UC(1)==0)
    E = 2*(rand(N,1)>fi)-1; % randomly defined E/I neurons
else
    % Generate grid structure excitatory/inhibitory
    %UC=[1 1 -1; -1 -1 -1; -1 -1 -1;] % unit cell
    %resSize=[11 11 1];
    sizeUC=size(UC); if(numel(sizeUC)<3); sizeUC(3)=1;end
    rpt=resSize./sizeUC;rpt=ceil(rpt);
    E=repmat(UC,rpt(1),rpt(2),rpt(3));E=E(1:resSize(1),1:resSize(2),1:resSize(3));
    E=E(:);
end

A = logical(sparse([],[],[],N,N));
N1 = resSize(1); N2 = resSize(2);N3 = resSize(3);
for n3=1:N3
    for n2 = 1:N2
        for n1 = 1:N1
            id = (n3-1)*N2*N1+(n2-1)*N1+n1;
            if(n3<N3)
                A(id,id+N1*N2) = 1;
            end
            if(n2<N2)
                A(id,id+N1) = 1;
            end
            if(n1<N1)
                A(id,id+1) = 1;
            end
            
            
            if(n3>1)
                A(id,id-N1*N2)= 1;
            end
            if(n2>1)
                A(id,id-N1)= 1;
            end
            if(n1>1)
                A(id,id-1) =1;
            end
            
        end
    end
end
G = 1.0*A;
G(E>0,E>0) = A(E>0,E>0)*W0(1,1); G(E>0,E<0) = A(E>0,E<0)*W0(1,2);
G(E<0,E>0) = A(E<0,E>0)*W0(2,1); G(E<0,E<0) = A(E<0,E<0)*W0(2,2);
[X,Xn,W] = find(G);
end

