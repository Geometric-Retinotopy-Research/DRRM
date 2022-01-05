function aic = calAIC(err)

err = double(err);


f = 0;
for i=1:size(err,1)
    
    mu = nanmean(err(i,:));
    sigma = nanstd(err(i,:));       
    for j =1:size(err,2)
        if g(err(i,j), mu, sigma) ==0
            continue
        end
        f = f + log(g(err(i,j), mu, sigma));
    end
end
k = size(err,1)*4;
aic = -2*f + 2*k;

aic = aic/size(err,1)/6;% we have six runs in each fitting

end

function y=g(x, mu, sigma)
y =  1/sigma/sqrt(2*pi)*exp(-0.5*(x-mu).^2/sigma^2);
if(isnan(y))
    y = 1e4;
end
end
 