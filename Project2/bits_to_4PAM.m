%This function takes a sequence of n bits and transforms it into 4PAM
%signal
function [S] = bits_to_4PAM(b)

%every time we check a pair of numbers
for k=1:2:length(b)
    if b(k)==0 && b(k+1)==0
        S(k) = 3;
    elseif(b(k)==0 && b(k+1)==1)
        S(k)=1;
    elseif(b(k)==1 && b(k+1)==1)    
        S(k)=-1;
    elseif(b(k)==1 && b(k+1)==0)    
        S(k)=-3;
    else
        disp('Error')
        return
    end
end
