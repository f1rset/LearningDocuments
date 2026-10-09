# 06. Диференціювання складної функції, матриці Якобі та Гессе (Multivariable Chain Rule, Jacobian & Hessian)

[⬅️ Попередня: 05. Функції багатьох змінних, частинні похідні та градієнт](05-multivariable-functions-limits-and-partial-derivatives.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Умовна оптимізація та множники Лагранжа ➡️](07-constrained-optimization-and-lagrange-multipliers.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Ланцюгове правило для багатьох змінних (Multivariable Chain Rule):** Похідні складних функцій та дерево залежностей.
2. **Матриця Якобі (Jacobian Matrix $J_f$):** Повна матриця перших похідних для векторних функцій $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$.
3. **Матриця Гессе (Hessian Matrix $H_f$):** Симетрична матриця других похідних, що описує кривизну поверхні.
4. **Багатовимірний ряд Тейлора (Multivariate Taylor Expansion):** Квадратична апроксимація функцій через градієнт та Гесіан.
5. **Класифікація екстремумів (Second Derivative Test):** Критерій Сильвестра, сідлові точки (Saddle points), локальні мінімуми та максимуми.

---

## 1. Ланцюгове правило (Multivariable Chain Rule)

### 1.1. Випадок однієї незалежної змінної: $z = f(x(t), y(t))$
Якщо $x = x(t)$ та $y = y(t)$ диференційовні за $t$, то повна похідна $\frac{dz}{dt}$ дорівнює:
$$\frac{dz}{dt} = \frac{\partial f}{\partial x}\frac{dx}{dt} + \frac{\partial f}{\partial y}\frac{dy}{dt} = \nabla f \cdot \mathbf{r}'(t)$$

### 1.2. Загальний випадок: $z = f(u(s, t), v(s, t))$
Кожна частинна похідна рахується як сума внесків за всіма проміжними шляхами:
$$\frac{\partial z}{\partial s} = \frac{\partial z}{\partial u}\frac{\partial u}{\partial s} + \frac{\partial z}{\partial v}\frac{\partial v}{\partial s}$$
$$\frac{\partial z}{\partial t} = \frac{\partial z}{\partial u}\frac{\partial u}{\partial t} + \frac{\partial z}{\partial v}\frac{\partial v}{\partial t}$$

```text
               z
             /   \
            u     v
           / \   / \
          s   t s   t
   ∂z/∂s = (∂z/∂u * ∂u/∂s) + (∂z/∂v * ∂v/∂s)
```

---

## 2. Матриця Якобі (Jacobian Matrix)

### 2.1. Означення для вектор-функції
Нехай $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$ перетворює вектор $\mathbf{x} = (x_1, \dots, x_n)^T$ у $\mathbf{f}(\mathbf{x}) = (f_1(\mathbf{x}), \dots, f_m(\mathbf{x}))^T$.
**Матриця Якобі (Jacobian Matrix $J_{\mathbf{f}} \in \mathbb{R}^{m \times n}$)** містить усі перші частинні похідні:

$$J_{\mathbf{f}}(\mathbf{x}) = \frac{\partial(f_1, \dots, f_m)}{\partial(x_1, \dots, x_n)} = \begin{pmatrix} 
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \dots & \frac{\partial f_1}{\partial x_n} \\ 
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \dots & \frac{\partial f_2}{\partial x_n} \\ 
\vdots & \vdots & \ddots & \vdots \\ 
\frac{\partial f_m}{\partial x_1} & \frac{\partial f_m}{\partial x_2} & \dots & \frac{\partial f_m}{\partial x_n} 
\end{pmatrix}$$

* **Матрична форма диференціала:** $d\mathbf{f} = J_{\mathbf{f}}(\mathbf{x})\,d\mathbf{x}$.
* **Матричне ланцюгове правило:** Для композиції $\mathbf{h}(\mathbf{x}) = \mathbf{g}(\mathbf{f}(\mathbf{x}))$:
  $$J_{\mathbf{h}}(\mathbf{x}) = J_{\mathbf{g}}(\mathbf{f}(\mathbf{x})) \cdot J_{\mathbf{f}}(\mathbf{x})$$

### 2.2. Якобіан (Jacobian Determinant)
Для відображення $n \to n$ визначник матриці Якобі $|J| = \det(J_{\mathbf{f}})$ називається **якобіаном**. Він показує коефіцієнт локальної зміни об'єму при заміні змінних:
$$dV_{\mathbf{y}} = |\det J_{\mathbf{f}}(\mathbf{x})|\,dV_{\mathbf{x}}$$

---

## 3. Матриця Гессе (Hessian Matrix)

### 3.1. Означення
Для скалярної функції $f: \mathbb{R}^n \to \mathbb{R}$ двічі неперервно диференційовної ($C^2$), **матриця Гессе (Hessian Matrix $H_f \in \mathbb{R}^{n \times n}$)** — це матриця других частинних похідних:

$$H_f(\mathbf{x}) = \nabla^2 f(\mathbf{x}) = \begin{pmatrix} 
\frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1 \partial x_2} & \dots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\ 
\frac{\partial^2 f}{\partial x_2 \partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \dots & \frac{\partial^2 f}{\partial x_2 \partial x_n} \\ 
\vdots & \vdots & \ddots & \vdots \\ 
\frac{\partial^2 f}{\partial x_n \partial x_1} & \frac{\partial^2 f}{\partial x_n \partial x_2} & \dots & \frac{\partial^2 f}{\partial x_n^2} 
\end{pmatrix}$$

* **Симетрія:** За теоремою Шварца матриця Гессе завжди **симетрична** ($H_f = H_f^T$).

---

## 4. Багатовимірна формула Тейлора (Multivariate Taylor Series)

Квадратичне наближення функції $f(\mathbf{x})$ в околі точки $\mathbf{x}_0$ з малим приростом $\mathbf{h} = \mathbf{x} - \mathbf{x}_0$:

$$f(\mathbf{x}_0 + \mathbf{h}) = f(\mathbf{x}_0) + \nabla f(\mathbf{x}_0)^T \mathbf{h} + \frac{1}{2} \mathbf{h}^T H_f(\mathbf{x}_0) \mathbf{h} + o(\|\mathbf{h}\|^2)$$

* **Лінійна частина:** $\nabla f(\mathbf{x}_0)^T \mathbf{h}$ (дотична площина).
* **Квадратична частина:** $\frac{1}{2} \mathbf{h}^T H_f(\mathbf{x}_0) \mathbf{h}$ (кривизна / чаша).

---

## 5. Дослідження на безумовний локальний екстремум (Extrema Analysis)

### 5.1. Необхідна умова екстремуму 1-го порядку (Fermat's Theorem)
У точці локального екстремуму $\mathbf{x}^*$:
$$\nabla f(\mathbf{x}^*) = \mathbf{0} \iff \frac{\partial f}{\partial x_1} = 0, \; \frac{\partial f}{\partial x_2} = 0, \; \dots, \; \frac{\partial f}{\partial x_n} = 0$$
Точки, де $\nabla f = \mathbf{0}$ або не існує, називаються **стаціонарними / критичними точками (Critical Points)**.

### 5.2. Достатня умова екстремуму для функції 2 змінних $f(x, y)$
Обчислюємо другий диференціал у стаціонарній точці $(x_0, y_0)$:
$$A = f_{xx}(x_0, y_0), \quad B = f_{xy}(x_0, y_0), \quad C = f_{yy}(x_0, y_0)$$
$$\det H = D = A C - B^2 = f_{xx} f_{yy} - (f_{xy})^2$$

1. **Якщо $D > 0$ та $A > 0$:** Точка $(x_0, y_0)$ є точкою **строгого локального мінімуму (Local Minimum)** (поверхня має форму параболічної чаші вгору).
2. **Якщо $D > 0$ та $A < 0$:** Точка $(x_0, y_0)$ є точкою **строгого локального максимуму (Local Maximum)** (чаша вниз).
3. **Якщо $D < 0$:** Точка $(x_0, y_0)$ є **сідловою точкою (Saddle Point)** (в одному напрямку мінімум, в іншому — максимум, екстремуму немає).
4. **Якщо $D = 0$:** Критерій не дає відповіді (потрібен аналіз вищих порядків).

```text
    Локальний мінімум (D > 0, A > 0)    Сідлова точка (D < 0)
              \     /                           ^ z
               \   /                            |   (Сідло)
                \_/                             +-------> y
             (x0, y0)                          /
                                            x /
```

### 5.3. Загальний випадок у $\mathbb{R}^n$ (Критерій Сильвестра / Власні значення)
Тип стаціонарної точки визначається додатною/від'ємною визначеністю матриці Гессе $H = H_f(\mathbf{x}^*)$:
* **Локальний мінімум:** $H \succ 0$ (додатно визначена: усі кутові мінори $\Delta_k > 0$ або всі власні значення $\lambda_i > 0$).
* **Локальний максимум:** $H \prec 0$ (від'ємно визначена: чергування знаків $(-1)^k \Delta_k > 0$ або всі $\lambda_i < 0$).
* **Сідлова точка:** $H$ невизначена (має як додатні, так і від'ємні власні значення $\lambda_i$).

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Неправильний порядок матричного множення в Chain Rule:** Матриці Якобі множаться як $J_{\mathbf{g}}(\mathbf{f}(\mathbf{x})) \cdot J_{\mathbf{f}}(\mathbf{x})$ (розмірності $m \times k$ та $k \times n \to m \times n$). Множення у зворотному порядку неможливе або хибне!
2. **Перевірка лише $A > 0$ без $D > 0$:** Якщо $f_{xx} > 0$, але $D < 0$, точка є сідловою, а НЕ мінімумом! Завжди спершу перевіряйте знак дискримінанта $D = AC - B^2$.
3. **Забута симетрія матриці Гессе:** $f_{xy} = f_{yx}$ для гладких функцій, тому матриця Гессе завжди симетрична з дійсними власними значеннями.

---

[⬅️ Попередня: 05. Функції багатьох змінних, частинні похідні та градієнт](05-multivariable-functions-limits-and-partial-derivatives.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Умовна оптимізація та множники Лагранжа ➡️](07-constrained-optimization-and-lagrange-multipliers.md) | [⚡ Cheat Sheet](cheat-sheet.md)
