% simple implementation of corr: this is to avoid the matlab corr function
function retval = corr (x, y) 
    c = cov (x, y);
    s = std (x)' * std (y);
    retval = c ./ s;
    retval = retval(1,2);

end

 
