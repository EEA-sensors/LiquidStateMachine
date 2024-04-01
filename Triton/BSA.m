%%  Ben Spiking Algorithm (encode signals into spike trains)
% ref :http://users.elis.ugent.be/~bschrauw/publicaties/ijcnn.pdf
% input : is matrix with time instance vector
%            input(i,k) = i'th signal at k'th time instant
% filter : is array of filter in time same for all(with same sampling)
% threshold : BSA algorithm threshold between errors
function [output,residual] = BSA(input,filter,threshold)
if(nargin<3)
    threshold = 0.955;
end
K = length(input(1,:));
filter = reshape(filter,[1 length(filter)]); % shape filter according to input
output = zeros(size(input));
for k = 1:K
    error1 = 0;
    error2 = 0;
    for j = 1 :length(filter)
        if (k+j-1 <= K)
            error1 = error1 + abs(input(:,k+j-1) - filter(j));
            error2 = error2 + abs(input(:,k+j-1));
        end
    end
    
    id = threshold <= (error2 - error1);
    output(:,k)=0;
    
    if(~isempty(id))
        output(id,k)=1;
        for j = 1:length(filter)
            if k+j-1 <= length(input)
                input(id,k+j-1) =input(id,k+j-1)- filter(j);
            end
        end
    end
end
residual = input;%(1:samples);
