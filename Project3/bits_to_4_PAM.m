function [X] = bits_to_4_PAM(b,A)
    %counter so that we don't have zeros
    counter=1;
    %create array of all possible outcomes
    Y = [-3*A, -1*A, A , 3*A];
    %creating space from the beginning for length(b)/2 since we have 4PAM
    X=zeros(1,length(b)/2);
    for i=1:2:length(b)
        if(b(i)==0 && b(i+1)==0)
            X(counter) = Y(1);
        elseif(b(i)==0 && b(i+1)==1)
            X(counter) = Y(2); 
        elseif(b(i)==1 && b(i+1)==1)
            X(counter) = Y(3);
        elseif(b(i)==1 && b(i+1)==0)
            X(counter) = Y(4);
        end
        counter=counter+1;
    end
end