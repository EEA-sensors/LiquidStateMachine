% create connection matrix from input to reservoir
% In : number of neurons in input layer
%     Number of neurons reservoir  (N)
%     (optional) input fanout in reservoir (k)
function C = applied_current_matrix(Nin,Nres,k)
if(nargin<3) 
    k = 4;
end
C = logical(zeros(Nres,Nin));
for i = 1:Nin
    id = randperm(Nres,k);
    C(id,i) = 1;
end
end



