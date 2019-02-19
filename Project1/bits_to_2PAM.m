%This function takes a sequence of n bits and transforms it into 2PAM
%signal
function [S] = bits_to_2PAM(b)
%Determine final length of S (Warning otherwise)
S = 1:length(b);

for n=1:length(b)
    if b(n)==0
        S(n) = +1;
    elseif b(n)==1
        S(n) = -1;
    else
        disp('Error: Not a binary was given');
        return;
    end
end
