# ⚡ Інженерна шпаргалка з математичного аналізу (Calculus & ODE Cheat Sheet)

[🏠 Головний зміст](index.md) | [📚 Повний курс лекцій](notes/01-limits-continuity-and-derivatives.md)

---

## 📌 TL;DR — Квінтесенція курсу (Core Highlights)

* **Похідна ($f'$):** Миттєва швидкість зміни / нахил дотичної. $f'(x) = \lim_{\Delta x \to 0} \frac{f(x+\Delta x)-f(x)}{\Delta x}$.
* **Інтеграл ($\int$):** Накопичення / орієнтована площа / обернена дія до диференціювання. $\int_a^b f(x)dx = F(b) - F(a)$.
* **Градієнт ($\nabla f$):** Вектор у напрямку найшвидшого зростання функції, ортогональний до ліній рівня. $\|\nabla f\| = \max D_{\mathbf{u}} f$.
* **Матриця Гессе ($H_f$):** Матриця кривизни поверхні. Додатно визначена ($H \succ 0$) $\implies$ локальний мінімум; від'ємно визначена ($H \prec 0$) $\implies$ локальний максимум; невизначена $\implies$ сідло.
* **Множники Лагранжа ($\nabla f = \lambda \nabla g$):** Дотик лінії рівня до межі обмеження.
* **Кратні інтеграли:** $dx\,dy = r\,dr\,d\theta$ (полярні), $dx\,dy\,dz = \rho^2 \sin\phi\,d\rho\,d\phi\,d\theta$ (сферичні).
* **Теореми векторного аналізу:** 
  * Грін: $\oint Pdx + Qdy = \iint (\partial_x Q - \partial_y P)dA$.
  * Гаусс-Остроградський: $\oiint \mathbf{F} \cdot d\mathbf{S} = \iiint (\nabla \cdot \mathbf{F}) dV$.
  * Стокс: $\oint \mathbf{F} \cdot d\mathbf{r} = \iint (\nabla \times \mathbf{F}) \cdot d\mathbf{S}$.
* **ДР 1-го порядку:** $y' + Py = Q \implies \mu = e^{\int P dx}, \, y = \frac{1}{\mu}\int \mu Q dx$.
* **ДР 2-го порядку:** $a y'' + b y' + c y = 0 \implies a\lambda^2 + b\lambda + c = 0$.

---

## 1. Таблиця похідних та інтегралів (Derivatives & Integrals Master Table)

| Функція $f(x)$ | Похідна $f'(x)$ | Первісна $\int f(x)\,dx$ |
| :--- | :--- | :--- |
| $x^n$ | $n x^{n-1}$ | $\frac{x^{n+1}}{n+1} + C \quad (n \neq -1)$ |
| $\frac{1}{x}$ | $-\frac{1}{x^2}$ | $\ln\|x\| + C$ |
| $e^x$ | $e^x$ | $e^x + C$ |
| $a^x$ | $a^x \ln a$ | $\frac{a^x}{\ln a} + C$ |
| $\sin x$ | $\cos x$ | $-\cos x + C$ |
| $\cos x$ | $-\sin x$ | $\sin x + C$ |
| $\tan x$ | $\frac{1}{\cos^2 x} = \sec^2 x$ | $-\ln\|\cos x\| + C$ |
| $\cot x$ | $-\frac{1}{\sin^2 x} = -\csc^2 x$ | $\ln\|\sin x\| + C$ |
| $\arcsin x$ | $\frac{1}{\sqrt{1 - x^2}}$ | $x \arcsin x + \sqrt{1 - x^2} + C$ |
| $\arctan x$ | $\frac{1}{1 + x^2}$ | $x \arctan x - \frac{1}{2}\ln(1 + x^2) + C$ |
| $\frac{1}{x^2 + a^2}$ | — | $\frac{1}{a} \arctan\left(\frac{x}{a}\right) + C$ |
| $\frac{1}{\sqrt{a^2 - x^2}}$ | — | $\arcsin\left(\frac{x}{a}\right) + C$ |
| $\frac{1}{x^2 - a^2}$ | — | $\frac{1}{2a} \ln\left\|\frac{x - a}{x + a}\right\| + C$ |
| $\frac{1}{\sqrt{x^2 \pm a^2}}$ | — | $\ln\left\|x + \sqrt{x^2 \pm a^2}\right\| + C$ |

---

## 2. Розклади Маклорена / Тейлора (Taylor Series Master Expansions)

| Функція | Розклад Маклорена ($x_0 = 0$) | Радіус збіжності $R$ |
| :--- | :--- | :--- |
| $\frac{1}{1 - x}$ | $1 + x + x^2 + x^3 + \dots = \sum_{n=0}^\infty x^n$ | $R = 1$ |
| $e^x$ | $1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \dots = \sum_{n=0}^\infty \frac{x^n}{n!}$ | $R = \infty$ |
| $\sin x$ | $x - \frac{x^3}{3!} + \frac{x^5}{5!} - \dots = \sum_{n=0}^\infty \frac{(-1)^n x^{2n+1}}{(2n+1)!}$ | $R = \infty$ |
| $\cos x$ | $1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \dots = \sum_{n=0}^\infty \frac{(-1)^n x^{2n}}{(2n)!}$ | $R = \infty$ |
| $\ln(1 + x)$ | $x - \frac{x^2}{2} + \frac{x^3}{3} - \dots = \sum_{n=1}^\infty \frac{(-1)^{n-1} x^n}{n}$ | $R = 1, \; (-1, 1]$ |
| $(1 + x)^\alpha$ | $1 + \alpha x + \frac{\alpha(\alpha-1)}{2!}x^2 + \dots = \sum_{n=0}^\infty \binom{\alpha}{n} x^n$ | $R = 1$ |

---

## 3. Матриця ознак збіжності числових рядів (Series Convergence Tests)

| Ознака | Формула / Умова | Збігається | Розбігається |
| :--- | :--- | :--- | :--- |
| **Необхідна** | $\lim_{n \to \infty} a_n$ | Може збігатися (якщо $0$) | $\lim a_n \neq 0$ |
| **$p$-Ряд** | $\sum \frac{1}{n^p}$ | $p > 1$ | $p \le 1$ |
| **Д'Аламбера** | $D = \lim \frac{a_{n+1}}{a_n}$ | $D < 1$ | $D > 1$ ($D=1$ — ?) |
| **Коші (корінь)** | $C = \lim \sqrt[n]{a_n}$ | $C < 1$ | $C > 1$ ($C=1$ — ?) |
| **Лейбніца** | $\sum (-1)^{n-1} b_n$ | $b_{n+1} \le b_n$ та $\lim b_n = 0$ | Не задовольняє |
| **Радіус $R$** | $R = \lim \left\|\frac{c_n}{c_{n+1}}\right\|$ | $\|x - x_0\| < R$ | $\|x - x_0\| > R$ |

---

## 4. Багатовимірний аналіз та оптимізація (Multivariable & Optimization)

### 4.1. Градієнт, похідна за напрямком та лінеаризація
$$\nabla f = \begin{pmatrix} \frac{\partial f}{\partial x_1}, \dots, \frac{\partial f}{\partial x_n} \end{pmatrix}^T, \quad D_{\mathbf{u}} f = \nabla f \cdot \mathbf{u} \quad (\|\mathbf{u}\| = 1)$$
$$f(\mathbf{x}_0 + \mathbf{h}) \approx f(\mathbf{x}_0) + \nabla f(\mathbf{x}_0)^T \mathbf{h} + \frac{1}{2}\mathbf{h}^T H_f(\mathbf{x}_0)\mathbf{h}$$

### 4.2. Класифікація стаціонарних точок ($\nabla f = \mathbf{0}$) у 2D
$$D = \det H = f_{xx} f_{yy} - (f_{xy})^2$$
* $D > 0$ та $f_{xx} > 0 \implies$ **Строгий локальний мінімум (Local Min)**.
* $D > 0$ та $f_{xx} < 0 \implies$ **Строгий локальний максимум (Local Max)**.
* $D < 0 \implies$ **Сідлова точка (Saddle Point)** (екстремуму немає).
* $D = 0 \implies$ Потрібен аналіз вищих порядків.

### 4.3. Умовна оптимізація Лагранжа (Lagrange Multipliers)
$$\min f(\mathbf{x}) \quad \text{s.t.} \quad g(\mathbf{x}) = 0 \implies \begin{cases} \nabla f = \lambda \nabla g \\ g(\mathbf{x}) = 0 \end{cases}$$

---

## 5. Заміни координат у кратних інтегралах (Coordinate Systems & Jacobians)

| Система | Зв'язок координат | Якобіан $|J|$ | Елемент об'єму/площі |
| :--- | :--- | :--- | :--- |
| **Декартові 2D** | $(x, y)$ | $1$ | $dA = dx\,dy$ |
| **Полярні 2D** | $x = r\cos\theta, \; y = r\sin\theta$ | $r$ | $dA = r\,dr\,d\theta$ |
| **Циліндричні 3D**| $x = r\cos\theta, \; y = r\sin\theta, \; z = z$ | $r$ | $dV = r\,dr\,d\theta\,dz$ |
| **Сферичні 3D** | $x = \rho\sin\phi\cos\theta, y = \rho\sin\phi\sin\theta, z = \rho\cos\phi$ | $\rho^2 \sin\phi$ | $dV = \rho^2 \sin\phi\,d\rho\,d\phi\,d\theta$ |

---

## 6. Інтегральні теореми векторного аналізу (Vector Calculus Theorems)

| Теорема | Формула | Сенс |
| :--- | :--- | :--- |
| **Ньютона-Лейбніца** | $\int_a^b F'(x)\,dx = F(b) - F(a)$ | 1D накопичення |
| **Потенціальне поле** | $\int_A^B \nabla \Phi \cdot d\mathbf{r} = \Phi(B) - \Phi(A)$ | Незалежність від шляху ($\operatorname{rot}\mathbf{F} = \mathbf{0}$) |
| **Гріна** (2D) | $\oint_{\partial D} Pdx + Qdy = \iint_D \left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right)dx\,dy$ | Контур $\to$ Площа |
| **Стокса** (3D) | $\oint_{\partial S} \mathbf{F} \cdot d\mathbf{r} = \iint_S (\operatorname{rot} \mathbf{F}) \cdot d\mathbf{S}$ | 3D Контур $\to$ Поверхня |
| **Остроградського-Гаусса** | $\oiint_{\partial V} \mathbf{F} \cdot d\mathbf{S} = \iiint_V (\operatorname{div} \mathbf{F})\,dV$ | Потік через оболонку = Джерела всередині |

---

## 7. Рецепти розв'язання диференціальних рівнянь (ODE Solving Recipes)

### 7.1. 1-й порядок
* **Відокремлювані:** $y' = f(x)g(y) \implies \int \frac{dy}{g(y)} = \int f(x)\,dx + C$.
* **Однорідні:** $y' = f(y/x) \implies u = y/x, \; y' = u'x + u$.
* **Лінійні:** $y' + P(x)y = Q(x) \implies \mu(x) = e^{\int P(x)dx} \implies y = \frac{1}{\mu(x)}\left(\int \mu(x)Q(x)dx + C\right)$.
* **Бернуллі:** $y' + Py = Q y^n \implies z = y^{1-n} \implies \frac{1}{1-n}z' + Pz = Q$.
* **Точні (в повних диф.):** $M dx + N dy = 0 \iff M_y = N_x \implies U(x, y) = C$.

### 7.2. 2-й порядок зі сталими коефіцієнтами: $a y'' + b y' + c y = 0$
Характеристичне рівняння: $a \lambda^2 + b \lambda + c = 0$:
1. $D > 0 \implies \lambda_1 \neq \lambda_2 \in \mathbb{R} \implies y = C_1 e^{\lambda_1 x} + C_2 e^{\lambda_2 x}$.
2. $D = 0 \implies \lambda_1 = \lambda_2 = \lambda \implies y = (C_1 + C_2 x) e^{\lambda x}$.
3. $D < 0 \implies \lambda = \alpha \pm i\beta \implies y = e^{\alpha x} (C_1 \cos(\beta x) + C_2 \sin(\beta x))$.

### 7.3. Класифікація фазових портретів $\dot{\mathbf{x}} = A\mathbf{x}$
* $\det A < 0 \implies$ **Сідлова точка (Saddle)** (завжди нестійка).
* $\det A > 0, D \ge 0 \implies$ **Вузол (Node)** (стійкий якщо $\operatorname{Tr} < 0$, нестійкий якщо $\operatorname{Tr} > 0$).
* $\det A > 0, D < 0, \operatorname{Tr} \neq 0 \implies$ **Фокус / Спіраль (Spiral)**.
* $\det A > 0, \operatorname{Tr} = 0 \implies$ **Центр (Center)** (еліпси, нейтральна стійкість).

---

[🏠 Головний зміст](index.md) | [📚 Повний курс лекцій](notes/01-limits-continuity-and-derivatives.md)

