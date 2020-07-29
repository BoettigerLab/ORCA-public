function a = HammingBound(n, d, q)

binomial = @(n,p) factorial(n)./(factorial(n-p).*factorial(p));

t = floor((d-1)/2);

k=0:t;

a =  q.^n./sum(binomial(n,k).*(q-1).^k);