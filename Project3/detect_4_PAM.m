function [PAMsymbols] = detect_4_PAM(data, A)
%we need to find which is point closer to our data input
X = [-3*A, -1*A, A , 3*A];

%initialize with zeros
PAMsymbols=zeros(1,length(data));

for i=1:length(data)
    %we find for each element of data the distance from each X
    ds1 = norm(X(1)-data(1,i));
    ds2 = norm(X(2)-data(1,i));
    ds3 = norm(X(3)-data(1,i));
    ds4 = norm(X(4)-data(1,i));
    
    %find the minimum value
    min_val = min([ds1,ds2,ds3,ds4]);
    %check every time which distance is the shortest
    if(ds1== min_val)
        PAMsymbols(1,i) = X(1);
    elseif(ds2 == min_val)
        PAMsymbols(1,i) = X(2);
    elseif(ds3 == min_val)
        PAMsymbols(1,i) = X(3);
    elseif(ds4 == min_val)
        PAMsymbols(1,i) = X(4);
    end
end

end