# 01. Границі, неперервність та диференційовність (Limits, Continuity & Single-Variable Derivatives)

[🏠 Головний зміст](index.md) | [Наступна тема: 02. Теореми аналізу та ряди Тейлора ➡️](02-mean-value-theorems-taylor-series-and-extrema.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Границя функції ($\varepsilon$-$\delta$ Definition of Limit):** Точна поведінка функції біля точки.
2. **Неперервність та розриви (Continuity & Classification of Discontinuities):** Усувні, 1-го та 2-го роду.
3. **Похідна (Derivative $f'(x)$):** Миттєва швидкість зміни та кутовий коефіцієнт дотичної.
4. **Правила диференціювання (Differentiation Rules):** Product rule, Quotient rule, Chain rule.
5. **Диференціал функції (Differential $df$):** Лінійна апроксимація функції.

---

## 1. Границі та неперервність (Limits & Continuity)

### 1.1. Означення границі за Коші ($\varepsilon$-$\delta$)
Нехай $f: D \to \mathbb{R}$. Число $L$ називається **границею функції $f(x)$ при $x \to x_0$** ($\lim_{x \to x_0} f(x) = L$), якщо:
$$\forall \varepsilon > 0 \quad \exists \delta > 0 \quad \forall x \in D: \quad 0 < |x - x_0| < \delta \implies |f(x) - L| < \varepsilon$$

> **💡 Інтуїція ("на пальцях"):**
> Наскільки б вузький горизонтальний коридор шириною $\pm \varepsilon$ навколо значення $L$ нам не задали, ми завжди можемо знайти настільки малий окіл $\pm \delta$ навколо $x_0$, що графік функції всередині цього околу не вийде за межі коридору.

### 1.2. Чудові границі (Important Standard Limits)
1. **Перша чудова границя (First Remarkable Limit):**
   $$\lim_{x \to 0} \frac{\sin x}{x} = 1, \quad \lim_{x \to 0} \frac{\tan x}{x} = 1, \quad \lim_{x \to 0} \frac{\arcsin x}{x} = 1$$
2. **Друга чудова границя (Число $e$ / Second Remarkable Limit):**
   $$\lim_{x \to \infty} \left(1 + \frac{1}{x}\right)^x = e, \quad \lim_{x \to 0} (1 + x)^{1/x} = e \approx 2.71828$$
3. **Логарифмічні та показникові границі:**
   $$\lim_{x \to 0} \frac{\ln(1+x)}{x} = 1, \quad \lim_{x \to 0} \frac{e^x - 1}{x} = 1, \quad \lim_{x \to 0} \frac{(1+x)^\alpha - 1}{x} = \alpha$$

### 1.3. Неперервність функції (Continuity)
Функція $f(x)$ є **неперервною в точці $x_0$**, якщо:
1. Функція визначена в точці $x_0$.
2. Границя $\lim_{x \to x_0} f(x)$ існує.
3. $\lim_{x \to x_0} f(x) = f(x_0) \iff \lim_{\Delta x \to 0} \Delta f = 0$.

* **Класифікація точок розриву (Discontinuities):**
  * **Усувний розрив (Removable discontinuity):** $\lim_{x \to x_0^-} f(x) = \lim_{x \to x_0^+} f(x) \neq f(x_0)$.
  * **Розрив 1-го роду / Стрибок (Jump discontinuity):** Односторонні границі скінченні, але не рівні: $\lim_{x \to x_0^-} f(x) \neq \lim_{x \to x_0^+} f(x)$.
  * **Розрив 2-го роду (Essential / Infinite discontinuity):** Хоча б одна з односторонніх границь дорівнює $\pm\infty$ або не існує.

---

## 2. Похідна та диференційовність (Derivatives)

### 2.1. Означення похідної
**Похідна функції (Derivative $f'(x)$ або $\frac{df}{dx}$)** в точці $x$ — це границя відношення приросту функції до приросту аргументу, коли приріст аргументу прямує до нуля:
$$f'(x) = \frac{df}{dx} = \lim_{\Delta x \to 0} \frac{\Delta f}{\Delta x} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x} = \lim_{h \to 0} \frac{f(x + h) - f(x)}{h}$$

* **Геометричний зміст:** $f'(x_0) = \tan\alpha = k$ — кутовий коефіцієнт (нахил) дотичної прямої до графіка $y = f(x)$ у точці $(x_0, f(x_0))$.
  * Рівняння дотичної (Tangent line): $y - f(x_0) = f'(x_0)(x - x_0)$.
  * Рівняння нормалі (Normal line): $y - f(x_0) = -\frac{1}{f'(x_0)}(x - x_0)$.
* **Фізичний зміст:** Миттєва швидкість зміни величини $v(t) = \frac{ds}{dt}$, прискорення $a(t) = \frac{dv}{dt} = \frac{d^2s}{dt^2}$, струм $i(t) = \frac{dq}{dt}$.

```text
               y ^
                 |                  / (Дотична з нахилом f'(x0))
                 |           .----*
                 |         /     (x0, f(x0))
                 |       /
               --+-----+-----------------> x
                 0    x0
```

### 2.2. Правила диференціювання (Differentiation Rules)
1. **Лінійність:** $(c f(x) + g(x))' = c f'(x) + g'(x)$.
2. **Добуток (Product Rule / Leibniz Rule):**
   $$(f \cdot g)' = f' g + f g'$$
3. **Частка (Quotient Rule):**
   $$\left( \frac{f}{g} \right)' = \frac{f' g - f g'}{g^2} \quad (g(x) \neq 0)$$
4. **Ланцюгове правило / Похідна складеної функції (Chain Rule):**
   $$\frac{d}{dx} f(g(x)) = f'(g(x)) \cdot g'(x) \iff \frac{df}{dx} = \frac{df}{dg} \cdot \frac{dg}{dx}$$
5. **Похідна оберненої функції (Inverse Function Rule):**
   $$(f^{-1})'(y_0) = \frac{1}{f'(x_0)}, \quad \text{де } y_0 = f(x_0), \ f'(x_0) \neq 0$$

---

## 3. Таблиця основних похідних (Table of Standard Derivatives)

| Функція $f(x)$ | Похідна $f'(x)$ | Область визначення |
|---|---|---|
| $c$ (константа) | $0$ | $\mathbb{R}$ |
| $x^n$ | $n x^{n-1}$ | $x > 0$ або $x \in \mathbb{R}$ для $n \in \mathbb{N}$ |
| $e^x$ | $e^x$ | $\mathbb{R}$ |
| $a^x$ | $a^x \ln a$ | $a > 0, a \neq 1$ |
| $\ln x$ | $\frac{1}{x}$ | $x > 0$ |
| $\log_a x$ | $\frac{1}{x \ln a}$ | $x > 0$ |
| $\sin x$ | $\cos x$ | $\mathbb{R}$ |
| $\cos x$ | $-\sin x$ | $\mathbb{R}$ |
| $\tan x$ | $\frac{1}{\cos^2 x} = 1 + \tan^2 x$ | $x \neq \frac{\pi}{2} + \pi k$ |
| $\cot x$ | $-\frac{1}{\sin^2 x}$ | $x \neq \pi k$ |
| $\arcsin x$ | $\frac{1}{\sqrt{1 - x^2}}$ | $|x| < 1$ |
| $\arccos x$ | $-\frac{1}{\sqrt{1 - x^2}}$ | $|x| < 1$ |
| $\arctan x$ | $\frac{1}{1 + x^2}$ | $\mathbb{R}$ |
| $\sinh x$ | $\cosh x$ | $\mathbb{R}$ |
| $\cosh x$ | $\sinh x$ | $\mathbb{R}$ |

---

## 4. Диференціал функції (Differential)

**Диференціал (Differential $df$)** — це головна лінійна частина приросту функції:
$$df = f'(x) \cdot \Delta x = f'(x) dx$$
* **Лінійна апроксимація (Linear Approximation / Tangent Approximation):**
  $$f(x_0 + \Delta x) \approx f(x_0) + df = f(x_0) + f'(x_0)\Delta x$$
  *(Наприклад: $\sqrt{1+x} \approx 1 + \frac{1}{2}x$, $(1+x)^\alpha \approx 1 + \alpha x$, $\sin x \approx x$, $e^x \approx 1 + x$ при $x \to 0$).*

---

## 📊 Зведена таблиця (Summary Table)

| Поняття | Математичний вираз | Опис |
|---|---|---|
| Границя | $\lim_{x \to x_0} f(x) = L$ | $\forall \varepsilon > 0\ \exists \delta > 0: 0 < \|x-x_0\| < \delta \implies \|f(x)-L\| < \varepsilon$ |
| Неперервність | $\lim_{x \to x_0} f(x) = f(x_0)$ | Немає стрибків та розривів |
| Похідна | $f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}$ | Кутовий коефіцієнт дотичної |
| Ланцюгове правило | $(f(g(x)))' = f'(g(x)) g'(x)$ | Диференціювання вкладених функцій |
| Диференціал | $df = f'(x) dx$ | Найкраща лінійна апроксимація |

---

[🏠 Головний зміст](index.md) | [Наступна тема: 02. Теореми аналізу та ряди Тейлора ➡️](02-mean-value-theorems-taylor-series-and-extrema.md) | [⚡ Cheat Sheet](cheat-sheet.md)
