function number = NChooseKWithBounds(n,k)
if k > n
    number = 0;
else
    number = nchoosek(n,k);
end
