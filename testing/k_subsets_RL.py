def algorithm_u(ns, m):
  def visit(n, a):
    print(a[-n:])
    
  def f(mu, nu, sigma, n, a):
    if mu == 2:
      visit(n, a)
    else:
      f(mu - 1, nu - 1, (mu + sigma) % 2, n, a)
    if nu == mu + 1:
      a[mu] = mu - 1
      visit(n, a)
      while a[nu] > 0:
        a[nu] = a[nu] - 1
        visit(n, a)
    elif nu > mu + 1:
      if (mu + sigma) % 2 == 1:
        a[nu - 1] = mu - 1
      else:
        a[mu] = mu - 1
      if (a[nu] + sigma) % 2 == 1:
        b(mu, nu - 1, 0, n, a)
      else:
        f(mu, nu - 1, 0, n, a)
      while a[nu] > 0:
        a[nu] = a[nu] - 1
        if (a[nu] + sigma) % 2 == 1:
          b(mu, nu - 1, 0, n, a)
        else:
          f(mu, nu - 1, 0, n, a)
  def b(mu, nu, sigma, n, a):
    if nu == mu + 1:
      while a[nu] < mu - 1:
        visit(n, a)
        a[nu] = a[nu] + 1
      visit(n, a)
      a[mu] = 0
    elif nu > mu + 1:
      if (a[nu] + sigma) % 2 == 1:
        f(mu, nu - 1, 0, n, a)
      else:
        b(mu, nu - 1, 0, n, a)
      while a[nu] < mu - 1:
        a[nu] = a[nu] + 1
        if (a[nu] + sigma) % 2 == 1:
          f(mu, nu - 1, 0, n, a)
        else:
          b(mu, nu - 1, 0, n, a)
      if (mu + sigma) % 2 == 1:
        a[nu - 1] = 0
      else:
        a[mu] = 0
    if mu == 2:
      visit(n, a)
    else:
      b(mu - 1, nu - 1, (mu + sigma) % 2, n, a)
  n = len(ns)
  a = [0] * (n + 1)
  for j in range(1, m + 1):
    a[n - m + j] = j - 1
  return f(m, n, 0, n, a)

algorithm_u([1, 2, 3, 4], 3)
