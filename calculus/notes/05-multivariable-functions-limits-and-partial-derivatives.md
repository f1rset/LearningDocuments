# 05. Функції багатьох змінних, частинні похідні та градієнт (Multivariable Functions, Partial Derivatives & Gradient)

[⬅️ Попередня: 04. Числові, степеневі та ряди Фур'є](04-numerical-power-and-fourier-series.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 06. Диференціювання складної функції, матриці Якобі та Гессе ➡️](06-multivariable-chain-rule-jacobian-and-hessian.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Функції кількох змінних ($f: \mathbb{R}^n \to \mathbb{R}$):** Область визначення, лінії/поверхні рівня (Level curves / Contour lines).
2. **Границя та неперервність у $\mathbb{R}^n$ (Multivariate Limits & Continuity):** Траєкторна залежність границі.
3. **Частинні похідні (Partial Derivatives $\frac{\partial f}{\partial x_i}$):** Швидкість зміни вздовж координатних осей.
4. **Теорема Шварца / Клеро (Schwarz's Theorem on Mixed Partials):** Рівність мішаних похідних $\frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x}$.
5. **Повний диференціал (Total Differential $df$):** Дотична гіперплощина та лінеаризація.
6. **Градієнт ($\nabla f$) та похідна за напрямком ($D_{\mathbf{u}} f$):** Вектор найшвидшого зростання, ортогональний до ліній рівня.

---

## 1. Функції багатьох змінних (Multivariable Functions)

### 1.1. Означення та геометрична інтерпретація
Функція $n$ змінних $f: D \to \mathbb{R}$ ($D \subseteq \mathbb{R}^n$) ставить у відповідність вектору $\mathbf{x} = (x_1, x_2, \dots, x_n)$ дійсне число $z = f(\mathbf{x})$.

* **Графік функції 2 змінних:** Поверхня в $\mathbb{R}^3$: $\{(x, y, z) \mid z = f(x, y), (x, y) \in D\}$.
* **Лінії рівня (Level Curves / Contours):** Множина точок на площині $xy$, де функція набуває сталого значення:
  $$\{(x, y) \in \mathbb{R}^2 \mid f(x, y) = C\}$$

```text
       z ^   Поверхня z = f(x, y)
         |       .-----.
         |      /       \
         |     |         |
         +-----|---------|-----> y
        /       \       /
     x /         '-----'
     
     Лінії рівня на площині xy:
         y ^     / C=3 \
           |   (  C=2   )
           |     \ C=1 /
           +-------------> x
```

### 1.2. Границі функцій багатьох змінних (Limits in $\mathbb{R}^n$)
$$\lim_{\mathbf{x} \to \mathbf{x}_0} f(\mathbf{x}) = L \iff \forall \varepsilon > 0 \; \exists \delta > 0 : \quad 0 < \|\mathbf{x} - \mathbf{x}_0\| < \delta \implies |f(\mathbf{x}) - L| < \varepsilon$$

* **Критична відмінність від $\mathbb{R}^1$:** У 1D ми наближаємося лише зліва або справа. У 2D і $\mathbb{R}^n$ існує **нескінченно багато шляхів** наближення (прямі $y = kx$, параболи $y = kx^2$, спіралі тощо).
* **Доведення неіснування границі:** Якщо при наближенні вздовж різних траєкторій виходять різні значення, границя **не існує**.
  * *Приклад:* $\lim_{(x,y)\to (0,0)} \frac{xy}{x^2 + y^2}$.
    * Покладемо $y = kx$: $\lim_{x\to 0} \frac{x(kx)}{x^2 + k^2x^2} = \frac{k}{1 + k^2}$. Границя залежить від кутового коефіцієнта $k$, тому границі в точці $(0,0)$ не існує!

---

## 2. Частинні похідні (Partial Derivatives)

### 2.1. Означення
**Частинна похідна $f$ за змінною $x$** у точці $(x_0, y_0)$ обчислюється як звичайна 1D-похідна при фіксованому $y = y_0$:
$$\frac{\partial f}{\partial x}(x_0, y_0) = \lim_{\Delta x \to 0} \frac{f(x_0 + \Delta x, y_0) - f(x_0, y_0)}{\Delta x}$$
Аналогічно за $y$ (фіксуючи $x = x_0$):
$$\frac{\partial f}{\partial y}(x_0, y_0) = \lim_{\Delta y \to 0} \frac{f(x_0, y_0 + \Delta y) - f(x_0, y_0)}{\Delta y}$$

* **Позначення:** $\frac{\partial f}{\partial x} = f'_x = f_x = \partial_x f$.

### 2.2. Мішані похідні та теорема Шварца (Schwarz's / Clairaut's Theorem)
Похідні вищих порядків позначаються:
$$\frac{\partial^2 f}{\partial x^2} = f_{xx}, \quad \frac{\partial^2 f}{\partial y^2} = f_{yy}, \quad \frac{\partial^2 f}{\partial x \partial y} = f_{yx} = \frac{\partial}{\partial x}\left(\frac{\partial f}{\partial y}\right)$$

> **👑 Теорема Шварца (Рівність мішаних похідних):**
> Якщо мішані частинні похідні $f_{xy}$ та $f_{yx}$ визначені та **неперервні** в околі точки $(x_0, y_0)$, то вони збігаються:
> $$\frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x}$$

---

## 3. Диференційовність та повний диференціал (Total Differential)

### 3.1. Означення диференційовності
Функція $f(x, y)$ називається **диференційовною** в точці $(x_0, y_0)$, якщо її повний приріст можна подати у вигляді:
$$\Delta f = f(x_0 + \Delta x, y_0 + \Delta y) - f(x_0, y_0) = \frac{\partial f}{\partial x}\Delta x + \frac{\partial f}{\partial y}\Delta y + \alpha \Delta x + \beta \Delta y$$
де $\alpha, \beta \to 0$ при $(\Delta x, \Delta y) \to (0,0)$.

* **Повний диференціал (Total Differential $df$):**
  $$df = \frac{\partial f}{\partial x} dx + \frac{\partial f}{\partial y} dy = \sum_{i=1}^n \frac{\partial f}{\partial x_i} dx_i$$

### 3.2. Дотична площина та нормаль (Tangent Plane & Normal)
* **Рівняння дотичної площини** до поверхні $z = f(x, y)$ у точці $(x_0, y_0, z_0)$:
  $$z - z_0 = f_x(x_0, y_0)(x - x_0) + f_y(x_0, y_0)(y - y_0)$$
* **Вектор нормалі до поверхні (Normal vector $\mathbf{n}$):**
  $$\mathbf{n} = \left( f_x(x_0, y_0), \; f_y(x_0, y_0), \; -1 \right)$$
* Для неявно заданої поверхні $F(x, y, z) = 0$:
  $$F_x(x_0, y_0, z_0)(x - x_0) + F_y(x_0, y_0, z_0)(y - y_0) + F_z(x_0, y_0, z_0)(z - z_0) = 0$$
  де нормаль $\mathbf{n} = \nabla F = (F_x, F_y, F_z)$.

---

## 4. Градієнт та похідна за напрямком (Gradient & Directional Derivative)

### 4.1. Вектор Градієнта (Gradient Vector $\nabla f$)
**Градієнт (Gradient)** функції $f(x_1, \dots, x_n)$ — це вектор-стовпець (або вектор-рядок) частинних похідних першого порядку:
$$\nabla f = \operatorname{grad} f = \begin{pmatrix} \frac{\partial f}{\partial x_1} \\ \frac{\partial f}{\partial x_2} \\ \vdots \\ \frac{\partial f}{\partial x_n} \end{pmatrix}$$

### 4.2. Похідна за напрямком (Directional Derivative $D_{\mathbf{u}} f$)
Швидкість зміни функції $f$ у напрямку одиничного вектора $\mathbf{u} = (u_1, u_2)$ ($\|\mathbf{u}\| = 1$):
$$D_{\mathbf{u}} f(\mathbf{x}) = \lim_{t \to 0} \frac{f(\mathbf{x} + t\mathbf{u}) - f(\mathbf{x})}{t} = \nabla f(\mathbf{x}) \cdot \mathbf{u} = \|\nabla f\| \cos\theta$$
де $\theta$ — кут між вектором градієнта $\nabla f$ та вектором напрямку $\mathbf{u}$.

### 4.3. Фундаментальні властивості градієнта
1. **Напрямок найшвидшого зростання (Steepest Ascent):**
   * Максимальне значення $D_{\mathbf{u}} f$ досягається при $\theta = 0$ ($\mathbf{u} \parallel \nabla f$) і дорівнює **нормі градієнта $\|\nabla f\|$**.
   * Градієнт вказує напрямок **найшвидшого локального зростання** функції.
2. **Напрямок найшвидшого спадання (Steepest Descent / Gradient Descent):**
   * Мінімальне значення досягається в напрямку $-\nabla f$ (основа методу градієнтного спуску в ML/Optimization).
3. **Ортогональність до ліній рівня:**
   * Вектор $\nabla f(\mathbf{x}_0)$ **перпендикулярний (ортогональний)** до лінії рівня $f(x, y) = C$, що проходить через точку $\mathbf{x}_0$.

```text
               y ^
                 |         \  Лінія рівня f(x, y) = C
                 |          \
                 |           * (x0, y0)
                 |          / \
                 |         /   \  ∇f (Перпендикулярний дотичній)
                 +--------/-----\-----> x
```

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Існування частинних похідних НЕ гарантує неперервність або диференційовність:**
   Функція може мати обидві похідні $f_x$ та $f_y$ у точці $(0,0)$, але бути навіть розривною в цій точці (потрібна умова неперервності самих похідних $f_x, f_y$).
2. **Забута нормалізація вектора при обчисленні $D_{\mathbf{v}} f$:**
   Формула $D_{\mathbf{v}} f = \nabla f \cdot \mathbf{v}$ працює **лише для одиничного вектора** $\|\mathbf{v}\| = 1$. Якщо вектор не одиничний, спершу треба нормалізувати: $\mathbf{u} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$.
3. **Плутанина між градієнтом і лінією рівня:** Градієнт напрямлений перпендикулярно до лінії рівня, а не вздовж неї (вздовж лінії рівня $D_{\mathbf{u}} f = 0$).

---

[⬅️ Попередня: 04. Числові, степеневі та ряди Фур'є](04-numerical-power-and-fourier-series.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 06. Диференціювання складної функції, матриці Якобі та Гессе ➡️](06-multivariable-chain-rule-jacobian-and-hessian.md) | [⚡ Cheat Sheet](cheat-sheet.md)
