# 08. Кратні інтеграли та заміни координат (Multiple Integrals & Coordinate Transformations)

[⬅️ Попередня: 07. Умовна оптимізація та множники Лагранжа](07-constrained-optimization-and-lagrange-multipliers.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 09. Векторний аналіз, криволінійні та поверхневі інтеграли ➡️](09-vector-calculus-and-line-surface-integrals.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Подвійний та потрійний інтеграли (Double & Triple Integrals):** Інтегрування за мірою площі $dA$ та об'єму $dV$.
2. **Теорема Фубіні (Fubini's Theorem):** Зведення кратного інтеграла до повторних 1D інтегралів та зміна порядку інтегрування.
3. **Заміна змінних у кратних інтегралах (Change of Variables & Jacobian):** $dx\,dy = |\det J|\,du\,dv$.
4. **Полярні координати (Polar Coordinates):** Якобіан $r$.
5. **Циліндричні та сферичні координати (Cylindrical & Spherical Coordinates):** Якобіани $r$ та $\rho^2 \sin\phi$.
6. **Застосування (Applications):** Маса, центр мас, моменти інерції, об'єми тіл.

---

## 1. Подвійний інтеграл Рімана (Double Integral)

### 1.1. Означення та геометричний зміст
Розглянемо замкнену обмежену область $D \subset \mathbb{R}^2$, розбиту на $n$ комірок площею $\Delta A_i$.

$$\iint_D f(x, y)\,dA = \lim_{\operatorname{diam} \to 0} \sum_{i=1}^n f(x_i^*, y_i^*)\,\Delta A_i$$

* **Геометричний зміст:** Якщо $f(x, y) \ge 0$, подвійний інтеграл дорівнює **об'єму циліндричного бруса** під поверхнею $z = f(x, y)$ над областю $D$.
* **Площа області $D$:** $S(D) = \iint_D 1\,dx\,dy$.

```text
       z ^          z = f(x, y)
         |         .--------.
         |        /|       /|
         |       / |      / |
         |      *--+-----*  |
         |      |  |     |  |
         +------|--+-----|--+-----> y
        /       | /      | /
     x /        '--------'  Область D на площині xy
```

### 1.2. Теорема Фубіні (Fubini's Theorem)
Якщо область $D$ описується нерівностями $a \le x \le b$ та $g_1(x) \le y \le g_2(x)$ (правильна в напрямку $Oy$), то подвійний інтеграл зводиться до **повторного (iterated integral)**:

$$\iint_D f(x, y)\,dx\,dy = \int_a^b \left( \int_{g_1(x)}^{g_2(x)} f(x, y)\,dy \right) dx$$

* **Зміна порядку інтегрування:** Якщо область правильно спроєктована на вісь $Oy$ ($c \le y \le d$, $h_1(y) \le x \le h_2(y)$):
  $$\iint_D f(x, y)\,dx\,dy = \int_c^d \left( \int_{h_1(y)}^{h_2(y)} f(x, y)\,dx \right) dy$$

---

## 2. Заміна змінних у подвійному інтегралі (Change of Variables)

Якщо відображення $x = x(u, v), y = y(u, v)$ є гладким та взаємно однозначним:
$$\iint_D f(x, y)\,dx\,dy = \iint_{D^*} f(x(u, v), y(u, v)) \cdot |J(u, v)|\,du\,dv$$

де **визначник Якобі (Якобіан)**:
$$J(u, v) = \frac{\partial(x, y)}{\partial(u, v)} = \det \begin{pmatrix} \frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\ \frac{\partial y}{\partial u} & \frac{\partial y}{\partial v} \end{pmatrix} = \frac{\partial x}{\partial u}\frac{\partial y}{\partial v} - \frac{\partial x}{\partial v}\frac{\partial y}{\partial u}$$

### 2.1. Полярні координати (Polar Coordinates)
$$x = r\cos\theta, \quad y = r\sin\theta, \quad r \ge 0, \; \theta \in [0, 2\pi)$$

* **Якобіан переходу:**
  $$J = \det \begin{pmatrix} \cos\theta & -r\sin\theta \\ \sin\theta & r\cos\theta \end{pmatrix} = r\cos^2\theta + r\sin^2\theta = r$$
  $$dx\,dy = r\,dr\,d\theta$$
* **Формула подвійного інтеграла в полярних координатах:**
  $$\iint_D f(x, y)\,dx\,dy = \iint_{D^*} f(r\cos\theta, r\sin\theta)\, r\,dr\,d\theta$$
* *Класичний приклад (Інтеграл Пуассона / Гауса):*
  $$I = \int_{-\infty}^\infty e^{-x^2}\,dx \implies I^2 = \int_{-\infty}^\infty \int_{-\infty}^\infty e^{-(x^2+y^2)}\,dx\,dy = \int_0^{2\pi} d\theta \int_0^\infty e^{-r^2} r\,dr = \pi \implies I = \sqrt{\pi}$$

---

## 3. Потрійний інтеграл (Triple Integral)

$$\iiint_V f(x, y, z)\,dV = \lim_{\operatorname{diam} \to 0} \sum_{i=1}^n f(x_i^*, y_i^*, z_i^*)\,\Delta V_i$$

### 3.1. Декартові координати
$$\iiint_V f(x, y, z)\,dx\,dy\,dz = \int_a^b dx \int_{y_1(x)}^{y_2(x)} dy \int_{z_1(x, y)}^{z_2(x, y)} f(x, y, z)\,dz$$

### 3.2. Циліндричні координати (Cylindrical Coordinates)
Зручні для тіл з осьовою або радіальною симетрією (циліндри, параболоїди, конуси):
$$x = r\cos\theta, \quad y = r\sin\theta, \quad z = z$$
$$|J| = r \implies dV = r\,dr\,d\theta\,dz$$

### 3.3. Сферичні координати (Spherical Coordinates)
Зручні для сферично-симетричних тіл (кулі, півсфери, конуси з вершиною в початку координат):
$$x = \rho \sin\phi \cos\theta, \quad y = \rho \sin\phi \sin\theta, \quad z = \rho \cos\phi$$
де $\rho \in [0, \infty)$ — радіус, $\phi \in [0, \pi]$ — зенітний кут (кут з віссю $Oz$), $\theta \in [0, 2\pi)$ — азимутальний кут на площині $xy$.

* **Якобіан сферичних координат:**
  $$|J| = \rho^2 \sin\phi \implies dV = \rho^2 \sin\phi\,d\rho\,d\phi\,d\theta$$
* *Приклад (Об'єм кулі радіуса $R$):*
  $$V = \int_0^{2\pi} d\theta \int_0^\pi \sin\phi\,d\phi \int_0^R \rho^2\,d\rho = (2\pi) \cdot (-\cos\pi + \cos 0) \cdot \left(\frac{R^3}{3}\right) = 2\pi \cdot 2 \cdot \frac{R^3}{3} = \frac{4}{3}\pi R^3$$

```text
               z ^        . P(x, y, z)
                 |       /|
                 |  ρ   / |
                 |  \ φ/  |
                 |   \/   |
                 +---|----+---------> y
                /    |   /
             x /     |  / r
                     v /
                      \  θ
```

---

## 4. Фізичні та геометричні застосування (Physical Applications)

Нехай $\rho(x, y, z)$ — густина розподілу маси в тілі $V$:
1. **Загальна маса (Total Mass):** $M = \iiint_V \rho(x, y, z)\,dV$.
2. **Центр мас (Center of Mass $(\bar{x}, \bar{y}, \bar{z})$):**
   $$\bar{x} = \frac{1}{M}\iiint_V x \rho\,dV, \quad \bar{y} = \frac{1}{M}\iiint_V y \rho\,dV, \quad \bar{z} = \frac{1}{M}\iiint_V z \rho\,dV$$
3. **Момент інерції відносно осі $Oz$ (Moment of Inertia):**
   $$I_z = \iiint_V (x^2 + y^2) \rho(x, y, z)\,dV$$

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Забутий множник Якобіана (Якобіан $r$ або $\rho^2 \sin\phi$):**
   Найпоширеніша помилка — писати $dx\,dy \to dr\,d\theta$ замість $r\,dr\,d\theta$, або $dx\,dy\,dz \to d\rho\,d\phi\,d\theta$ замість $\rho^2 \sin\phi\,d\rho\,d\phi\,d\theta$.
2. **Плутанина між кутами $\phi$ та $\theta$ у сферичних координатах:**
   Зенітний кут $\phi$ змінюється від $0$ до $\pi$, тоді як азимутальний $\theta$ — від $0$ до $2\pi$.
3. **Неправильний порядок меж у повторному інтегралі:** Межі зовнішнього інтеграла **завжди** повинні бути сталими числами, а не функціями від внутрішніх змінних!

---

[⬅️ Попередня: 07. Умовна оптимізація та множники Лагранжа](07-constrained-optimization-and-lagrange-multipliers.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 09. Векторний аналіз, криволінійні та поверхневі інтеграли ➡️](09-vector-calculus-and-line-surface-integrals.md) | [⚡ Cheat Sheet](cheat-sheet.md)
