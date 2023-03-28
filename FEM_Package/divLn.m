function divLn(mark,n,title)
% Dividing line.
if nargin == 2
    disp(repmat(mark,[1,n]));
elseif nargin >= 3
    disp( [repmat(mark,[1,floor((n-2-length(title))/2)]), ' ', title, ' ', repmat(mark,[1,ceil((n-2-length(title))/2)])] );
end
end

