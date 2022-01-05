function index = getIndex(start, last, M)
index = zeros(1,M);
index(M) = last;
for i = 1:M-1
    if i==1
        index(i) = floor((last-1-start)/(M-1)) + start;
    else
       index(i) = floor((last-1-start)/(M-1)) + index(i-1);
    end
end
end
