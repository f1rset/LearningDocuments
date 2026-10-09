# 03. Невизначений та визначений інтеграли (Indefinite & Definite Integrals)

[⬅️ Попередня: 02. Теореми аналізу та ряди Тейлора](02-mean-value-theorems-taylor-series-and-extrema.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 04. Числові, степеневі та ряди Фур'є ➡️](04-numerical-power-and-fourier-series.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Первісна та невизначений інтеграл (Antiderivative & Indefinite Integral):** Операція, обернена до взяття похідної.
2. **Методи інтегрування (Integration Techniques):** Заміна змінних (Substitution), інтегрування частинами (Integration by parts), інтегрування раціональних дробів (Partial fraction decomposition).
3. **Визначений інтеграл Рімана (Definite Riemann Integral):** Границя інтегральних сум, орієнтована площа під кривою.
4. **Основна теорема аналізу (Fundamental Theorem of Calculus):** Формула Ньютона-Лейбніца.
5. **Невласні інтеграли (Improper Integrals):** Нескінченні межі та необмежені підінтегральні функції.

---

## 1. Первісна та невизначений інтеграл (Antiderivatives)

### 1.1. Означення
Функція $F(x)$ називається **первісною (Antiderivative)** для $f(x)$ на інтервалі $(a, b)$, якщо:
$$\forall x \in (a, b): \quad F'(x) = f(x)$$

Сукупність усіх первісних для функції $f(x)$ називається **невизначеним інтегралом (Indefinite Integral)**:
$$\int f(x)\,dx = F(x) + C, \quad C \in \mathbb{R}$$

> **💡 Інтуїція ("на пальцях"):**
> Диференціювання знаходить швидкість руху, знаючи координату. Інтегрування робить зворотне — знаючи швидкість у кожен момент часу $v(t)$, відновлює пройдену траєкторію $s(t)$ з точністю до початкової точки відліку $C$.

### 1.2. Таблиця базових інтегралів (Standard Integrals Table)

| Функція $f(x)$ | Невизначений інтеграл $\int f(x)\,dx$ | Умови |
| :--- | :--- | :--- |
| $x^n$ | $\frac{x^{n+1}}{n+1} + C$ | $n \neq -1$ |
| $\frac{1}{x}$ | $\ln\|x\| + C$ | $x \neq 0$ |
| $e^x$ | $e^x + C$ | — |
| $a^x$ | $\frac{a^x}{\ln a} + C$ | $a > 0, a \neq 1$ |
| $\sin x$ | $-\cos x + C$ | — |
| $\cos x$ | $\sin x + C$ | — |
| $\frac{1}{\cos^2 x} = \sec^2 x$ | $\tan x + C$ | $x \neq \frac{\pi}{2} + \pi k$ |
| $\frac{1}{\sin^2 x} = \csc^2 x$ | $-\cot x + C$ | $x \neq \pi k$ |
| $\frac{1}{x^2 + a^2}$ | $\frac{1}{a} \arctan\left(\frac{x}{a}\right) + C$ | $a \neq 0$ |
| $\frac{1}{\sqrt{a^2 - x^2}}$ | $\arcsin\left(\frac{x}{a}\right) + C$ | $\|x\| < a$ |
| $\frac{1}{x^2 - a^2}$ | $\frac{1}{2a}\ln\left\|\frac{x - a}{x + a}\right\| + C$ | $\|x\| \neq a$ ("високий логарифм") |
| $\frac{1}{\sqrt{x^2 \pm a^2}}$ | $\ln\left\|x + \sqrt{x^2 \pm a^2}\right\| + C$ | ("довгий логарифм") |

---

## 2. Основні методи інтегрування (Integration Techniques)

### 2.1. Метод заміни змінної (Substitution Rule / $u$-Substitution)
Якщо $x = \varphi(t)$ — диференційовна монотонна функція, то:
$$\int f(x)\,dx = \int f(\varphi(t))\,\varphi'(t)\,dt$$

* **Підведення під знак диференціала (Integration by inspection):**
  $$\int f(g(x))\,g'(x)\,dx = \int f(u)\,du, \quad \text{де } u = g(x), \, du = g'(x)\,dx$$
  * *Приклад:* $\int 2x e^{x^2}\,dx = \int e^{x^2}\,d(x^2) = e^{x^2} + C$.

### 2.2. Інтегрування частинами (Integration by Parts)
Випливає з правила похідної добутку $(uv)' = u'v + uv'$:
$$\int u\,dv = u v - \int v\,du$$

* **Правило вибору $u$ (Мнемоніка LIATE / LATE):**
  1. **L** — Logarithmic ($\ln x, \log_a x$)
  2. **I** — Inverse trigonometric ($\arcsin x, \arctan x$)
  3. **A** — Algebraic / Polynomials ($x^n, 3x^2 + 1$)
  4. **T** — Trigonometric ($\sin x, \cos x$)
  5. **E** — Exponential ($e^x, 2^x$)

  *Приклад:* $\int x \cos x\,dx$. Обираємо $u = x \implies du = dx$, $dv = \cos x\,dx \implies v = \sin x$.
  $$\int x \cos x\,dx = x \sin x - \int \sin x\,dx = x \sin x + \cos x + C$$

### 2.3. Інтегрування раціональних дробів (Partial Fractions)
Будь-який правильний раціональний дріб $\frac{P(x)}{Q(x)}$ ($\deg P < \deg Q$) розкладається на суму елементарних дробів:
$$\frac{P(x)}{(x-a)^k(x^2+px+q)^m} = \frac{A_1}{x-a} + \dots + \frac{A_k}{(x-a)^k} + \frac{M_1 x + N_1}{x^2+px+q} + \dots$$

---

## 3. Визначений інтеграл Рімана (Definite Riemann Integral)

### 3.1. Інтегральна сума та означення інтеграла
Розглянемо розбиття відрізка $[a, b]$ на $n$ частин: $a = x_0 < x_1 < \dots < x_n = b$.
Нехай $\Delta x_i = x_i - x_{i-1}$, $\xi_i \in [x_{i-1}, x_i]$, діаметр розбиття $\lambda = \max \Delta x_i$.

**Інтеграл Рімана** — це границя інтегральної суми:
$$\int_a^b f(x)\,dx = \lim_{\lambda \to 0} \sum_{i=1}^n f(\xi_i)\,\Delta x_i$$

```text
    y ^
      |           f(x)
      |          .---*---.
      |         /|   |   |\
      |        / |   |   | \
      |       /  |   |   |  \
      |      |   |   |   |   |
      +------|---|---|---|---|--------> x
      0      a  x1  xi  xn-1 b
                 |<-Δxi->|
```

* **Геометричний зміст:** Орієнтована (знакова) площа фігури, обмеженої графіком $y = f(x)$, віссю $Ox$ та прямими $x = a, x = b$.

### 3.2. Основна теорема аналізу (Fundamental Theorem of Calculus)
1. **Перша частина (Диференціювання інтеграла зі змінною верхньою межею):**
   $$\frac{d}{dx}\left(\int_a^x f(t)\,dt\right) = f(x)$$
2. **Друга частина — Формула Ньютона-Лейбніца (Newton-Leibniz Formula):**
   Якщо $f(x)$ неперервна на $[a, b]$, а $F(x)$ — її первісна ($F' = f$):
   $$\int_a^b f(x)\,dx = F(b) - F(a) = \left. F(x) \right|_a^b$$

### 3.3. Властивості визначеного інтеграла
1. **Лінійність (Linearity):** $\int_a^b (\alpha f(x) + \beta g(x))\,dx = \alpha \int_a^b f(x)\,dx + \beta \int_a^b g(x)\,dx$.
2. **Адитивність по відрізку (Additivity):** $\int_a^b f(x)\,dx = \int_a^c f(x)\,dx + \int_c^b f(x)\,dx$.
3. **Теорема про середнє значення (Mean Value Theorem for Integrals):**
   $$\exists c \in [a, b]: \quad \int_a^b f(x)\,dx = f(c)(b - a) \iff f(c) = \frac{1}{b-a}\int_a^b f(x)\,dx$$
4. **Заміна меж при інтегруванні за формулою Ньютона-Лейбніца:**
   $$\int_a^b f(\varphi(t))\varphi'(t)\,dt = \int_{\varphi(a)}^{\varphi(b)} f(u)\,du$$
5. **Інтегрування частинами у визначеному інтегралі:**
   $$\int_a^b u\,dv = \left. u v \right|_a^b - \int_a^b v\,du$$

---

## 4. Застосування визначеного інтеграла (Geometric & Physical Applications)

1. **Площа криволінійної трапеції (Area under curve):**
   $$S = \int_a^b |f(x) - g(x)|\,dx$$
2. **Довжина дуги кривої (Arc Length):**
   $$L = \int_a^b \sqrt{1 + (f'(x))^2}\,dx, \quad L = \int_{t_1}^{t_2} \sqrt{(x'(t))^2 + (y'(t))^2}\,dt$$
3. **Об'єм тіла обертання навколо осі $Ox$ (Volume of Revolution):**
   $$V_x = \pi \int_a^b (f(x))^2\,dx$$
4. **Площа поверхні обертання навколо осі $Ox$ (Surface Area of Revolution):**
   $$Q_x = 2\pi \int_a^b |f(x)|\sqrt{1 + (f'(x))^2}\,dx$$
5. **Робота сили $F(x)$ (Work Done):** $W = \int_{x_1}^{x_2} F(x)\,dx$.

---

## 5. Невласні інтеграли (Improper Integrals)

### 5.1. 1-го роду: Нескінченний проміжок інтегрування (Infinite interval)
$$\int_a^\infty f(x)\,dx = \lim_{B \to \infty} \int_a^B f(x)\,dx$$
Якщо границя скінченна — інтеграл **збігається (converges)**, якщо $\pm\infty$ або не існує — **розбігається (diverges)**.

* **Еталонний інтеграл:**
  $$\int_1^\infty \frac{dx}{x^p} = \begin{cases} \text{збігається}, & p > 1 \\ \text{розбігається}, & p \le 1 \end{cases}$$

### 5.2. 2-го роду: Необмежена функція (Unbounded integrand / Singular point)
Нехай $f(x)$ прямує до $\infty$ при $x \to b^-$:
$$\int_a^b f(x)\,dx = \lim_{\varepsilon \to 0^+} \int_a^{b - \varepsilon} f(x)\,dx$$

* **Еталонний інтеграл:**
  $$\int_a^b \frac{dx}{(b - x)^p} = \begin{cases} \text{збігається}, & p < 1 \\ \text{розбігається}, & p \ge 1 \end{cases}$$

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Забута стала $C$ у невизначеному інтегралі:** $\int f(x)\,dx = F(x)$ замість $F(x) + C$.
2. **Неперерахунок меж інтегрування при заміні змінних:** Якщо робиться заміна $u = g(x)$, межі обов'язково змінюються з $[a, b]$ на $[g(a), g(b)]$.
3. **Застосування формули Ньютона-Лейбніца до розривної функції:**
   $\int_{-1}^1 \frac{1}{x^2}\,dx \neq \left. -\frac{1}{x} \right|_{-1}^1 = -1 - 1 = -2$ ❌ (помилка, оскільки функція $1/x^2 > 0$ всюди, а в $x=0$ нескінченний розрив; правильний інтеграл $\int_{-1}^1 \frac{1}{x^2} dx = +\infty$ розбігається!).

---

[⬅️ Попередня: 02. Теореми аналізу та ряди Тейлора](02-mean-value-theorems-taylor-series-and-extrema.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 04. Числові, степеневі та ряди Фур'є ➡️](04-numerical-power-and-fourier-series.md) | [⚡ Cheat Sheet](cheat-sheet.md)

