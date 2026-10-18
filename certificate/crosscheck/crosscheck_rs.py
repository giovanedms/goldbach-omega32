"""Cross-check: sanidade numerica dos enunciados herdados (ate 1e8) + dados para fechar Sum log p/p < log x em x >= 1e8
via theta(t) < t (Buthe 2018, t <= 1e19) e Axler 2203.05917 Thm 1.6 (x >= 1.757e12).
Checa nos pontos criticos (x = p e x -> p^-): RS Thm 5 lado inferior (x>1) e superior (x>=286), Sum log p/p < log x."""
import math
import numpy as np

B = 0.2614972128476428
n = 10 ** 8
s = np.ones(n + 1, dtype=bool); s[:2] = False
for i in range(2, int(n ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
p = np.flatnonzero(s).astype(float)
S1 = np.cumsum(1 / p)                  # sum_{q<=p}
lp = np.log(p); llp = np.log(lp)
# x = p (sum includes p): upper bound worst at x=p ; lower bound worst at x -> next prime^- (sum incl p, x -> p_next)
up = S1 - (llp + B + 1 / (2 * lp ** 2))
m = p >= 286
print("RS Thm5 superior (x=p>=286): max(sum - rhs) =", up[m].max(), "(precisa <0)")
xn = p[1:]; lxn = np.log(xn)
low = S1[:-1] - (np.log(lxn) + B - 1 / (2 * lxn ** 2))
print("RS Thm5 inferior (x->p_next^-): min(sum - rhs) =", low.min(), "(precisa >0)")
S2 = np.cumsum(lp / p)
print("Sum log p/p - log x (x=p): max =", (S2 - lp).max(), "(precisa <0)")
th = lp.sum()
print(f"x0=1e8: Sum_(p<=x0) log p/p = {S2[-1]:.10f}, log x0 = {math.log(n):.10f}, theta(x0)/x0 = {th / n:.10f}")
print("  => para x>=x0 (theta(t)<t): Sum < Sum(x0) + 1 - theta(x0)/x0 + log(x/x0) = log x +",
      S2[-1] + 1 - th / n - math.log(n))
