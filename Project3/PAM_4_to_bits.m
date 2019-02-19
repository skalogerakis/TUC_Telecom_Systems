function [revBits] = PAM_4_to_bits(X,A)
%the opposite from bits_to_4_PAM
%transforms PAM_4 decoding to bits sequence

Y = [-3*A, -1*A, A , 3*A];

counter=1;
%again create default space
revBits = zeros(1,2*length(X));
for i=1:length(X)
    if(X(i) == Y(1))
        revBits(counter) =0;
        revBits(counter+1) = 0;
    elseif(X(i) == Y(2))
        revBits(counter) =0;
        revBits(counter+1) = 1;
    elseif(X(i) == Y(3))
        revBits(counter) =1;
        revBits(counter+1) = 1;
    elseif(X(i) == Y(4))
        revBits(counter) =1;
        revBits(counter+1) = 0;
    end
    %the counter increases by 2 for 
    counter = counter+2;
end

end