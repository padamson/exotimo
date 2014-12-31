for n = [2 4 9]
  disp(sprintf('factorial(%g)!=%g', n, factorial(n)));
end

for m = [10 15 20]
  for n = [3 5 7]
    disp(sprintf('binomial(%g,%g)=%g', m, n, binomial(m,n)));
  end
end

