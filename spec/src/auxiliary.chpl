proc factorial(n: int): int{
  if (n <= 1) then return 1;
  return n*factorial(n-1);
}
proc binomial(m: int, n: int): int {
  return factorial(m)/(factorial(n)*factorial(m-n));
}
