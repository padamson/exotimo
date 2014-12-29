proc factorial(n: int): int{
  if (n <= 1) then return 1;
  return n*factorial(n-1);
}
