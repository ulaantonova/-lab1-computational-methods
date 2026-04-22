# Модель: Метод Ньютона (5 семестр)
# Автор: Антонова Юлія, група АІ-231
# Функція та її похідна
def f(x):
return 2*x**3 - 9*x**2 - 60*x + 1
def df(x):
return 6*x**2 - 18*x - 60
# --- Бісекція ---
def bisection(func, a, b, eps=1e-4, max_iter=1000):
fa, fb = func(a), func(b)
if fa * fb > 0:
raise ValueError("На відрізку немає зміни знаку")
it = 0
while (b - a) / 2 > eps and it < max_iter:
c = (a + b) / 2
fc = func(c)
if fc == 0:
return c, it
if fa * fc < 0:
b, fb = c, fc
else:
a, fa = c, fc
it += 1
return (a + b) / 2, it
# --- Ньютон ---
def newton(func, dfunc, x0, eps=1e-4, max_iter=1000):
x = x0
for it in range(1, max_iter+1):
fx, dfx = func(x), dfunc(x)
if abs(dfx) < 1e-12:
return x, it, False
dx = -fx/dfx
x += dx
if abs(dx) < eps:
return x, it, True
return x, max_iter, False
x_vals = np.linspace(-10, 10, 1000)
y_vals = f(x_vals)
intervals = []
for i in range(len(x_vals)-1):
if y_vals[i] * y_vals[i+1] < 0:
intervals.append((x_vals[i], x_vals[i+1]))
results = []
for (a, b) in intervals:
# Бісекція
root_bis, it_bis = bisection(f, a, b)
# Ньютон
root_new, it_new, conv = newton(f, df, (a+b)/2)
results.append({
"interval": f"[{a:.3f}, {b:.3f}]",
"bis_root": round(root_bis, 6),
"bis_iters": it_bis,
"newt_root": round(root_new, 6),
"newt_iters": it_new,
"newt_conv": conv
})
df = pd.DataFrame(results)
print(df)
# --- Побудова графіка ---
plt.axhline(0, color='black')
plt.plot(x_vals, y_vals, label="f(x)")
for r in df["bis_root"]:
plt.plot(r, f(r), 'ro', label="Корінь (бісекція)")
plt.legend()
plt.grid(True)
plt.show()
# --- Запис у файл ---
with open("results_roots.txt", "w", encoding="utf-8") as fobj:
fobj.write("Розв'язок рівняння 2x^3 - 9x^2 - 60x + 1 = 0\n\n")
fobj.write(df.to_string(index=False))
print("Результати збережено у файл results_roots.txt")
