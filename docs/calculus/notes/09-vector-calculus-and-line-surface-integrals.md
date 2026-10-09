# 09. Векторний аналіз, криволінійні та поверхневі інтеграли (Vector Calculus, Line & Surface Integrals)

[⬅️ Попередня: 08. Кратні інтеграли та заміни змінних](08-multiple-integrals-and-coordinate-transformations.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Диференціальні рівняння та системи ➡️](10-differential-equations-ode-and-systems.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Векторні поля та диференціальні оператори (Vector Fields & Differential Operators):** Градієнт ($\nabla f$), Дивергенція ($\nabla \cdot \mathbf{F}$), Ротор ($\nabla \times \mathbf{F}$), Лапласіан ($\Delta f$).
2. **Криволінійні інтеграли 1-го та 2-го роду (Line Integrals of 1st & 2nd Kind):** Довжина/маса кривої та робота силового поля ($\int_C \mathbf{F} \cdot d\mathbf{r}$).
3. **Потенціальні поля (Conservative Vector Fields):** Незалежність інтеграла від шляху інтегрування та потенціал поля ($\mathbf{F} = \nabla \Phi$).
4. **Теорема Гріна (Green's Theorem):** Зв'язок циркуляції по контуру з подвійним інтегралом на площині.
5. **Поверхневі інтеграли та потік векторного поля (Surface Integrals & Flux):** $\iint_S \mathbf{F} \cdot d\mathbf{S}$.
6. **Теорема Остроградського-Гаусса (Divergence Theorem):** Потік через замкнену поверхню дорівнює інтегралу від дивергенції по об'єму.
7. **Теорема Стокса (Stokes' Theorem):** Циркуляція по замкненому просторовому контуру дорівнює потоку ротора через натягнуту поверхню.

---

## 1. Диференціальні оператори векторного аналізу

Нехай $f(x, y, z)$ — скалярне поле, $\mathbf{F}(x, y, z) = P(x, y, z)\mathbf{i} + Q(x, y, z)\mathbf{j} + R(x, y, z)\mathbf{k}$ — векторне поле.
Оператор Гамільтона (набла): $\nabla = \mathbf{i}\frac{\partial}{\partial x} + \mathbf{j}\frac{\partial}{\partial y} + \mathbf{k}\frac{\partial}{\partial z}$.

| Оператор | Математичний запис | Фізичний зміст |
| :--- | :--- | :--- |
| **Градієнт (Gradient)** | $\nabla f = \left(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z}\right)$ | Напрямок найшвидшого зростання величини |
| **Дивергенція (Divergence)** | $\operatorname{div} \mathbf{F} = \nabla \cdot \mathbf{F} = \frac{\partial P}{\partial x} + \frac{\partial Q}{\partial y} + \frac{\partial R}{\partial z}$ | Густина джерел (+) або стоків (-) поля |
| **Ротор (Curl / Rotor)** | $\operatorname{rot} \mathbf{F} = \nabla \times \mathbf{F} = \begin{vmatrix} \mathbf{i} & \mathbf{j} & \mathbf{k} \\ \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z} \\ P & Q & R \end{vmatrix}$ | Локальна вихрова структура / закрученість поля |
| **Лапласіан (Laplacian)** | $\Delta f = \nabla^2 f = \operatorname{div}(\nabla f) = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} + \frac{\partial^2 f}{\partial z^2}$ | Нерівномірність розподілу поля відносно середнього |

* **Фундаментальні тотожності (Fundamental Identities):**
  1. $\operatorname{rot}(\operatorname{grad} f) = \nabla \times (\nabla f) = \mathbf{0}$ (Градієнтне поле завжди безвихрове).
  2. $\operatorname{div}(\operatorname{rot} \mathbf{F}) = \nabla \cdot (\nabla \times \mathbf{F}) = 0$ (Вихрове поле не має джерел і стоків).

---

## 2. Криволінійні інтеграли (Line Integrals)

### 2.1. Інтеграл 1-го роду (за довжиною дуги $ds$)
Інтегрування скалярної функції $f(x, y, z)$ вздовж гладкої кривої $C$, параметризованої $\mathbf{r}(t) = (x(t), y(t), z(t)), t \in [a, b]$:
$$\int_C f(x, y, z)\,ds = \int_a^b f(x(t), y(t), z(t)) \cdot \|\mathbf{r}'(t)\|\,dt$$
де елемент дуги $ds = \sqrt{(x')^2 + (y')^2 + (z')^2}\,dt$.
* *Застосування:* Маса кривої $M = \int_C \rho(x, y, z)\,ds$, довжина дуги $L = \int_C 1\,ds$.

### 2.2. Інтеграл 2-го роду (за координатами / робота векторного поля)
Інтегрування векторного поля вздовж орієнтованої кривої $C$:
$$\int_C \mathbf{F} \cdot d\mathbf{r} = \int_C P\,dx + Q\,dy + R\,dz = \int_a^b \left( P x'(t) + Q y'(t) + R z'(t) \right) dt$$
* *Фізичний зміст:* **Робота** силового поля $\mathbf{F}$ при переміщенні матеріальної точки вздовж траєкторії $C$.
* *Властивість:* Зміна орієнтації кривої змінює знак: $\int_{C^-} \mathbf{F} \cdot d\mathbf{r} = -\int_{C^+} \mathbf{F} \cdot d\mathbf{r}$.

```text
               y ^
                 |         . B (x(b), y(b))
                 |        /   ^ F (Векторне поле)
                 |       /   /
                 |      /   /  dr (Дотичний вектор)
                 |   A .---'
                 +---------------------> x
                 Робота: W = ∫_C F · dr
```

---

## 3. Потенціальні векторні поля (Conservative Fields)

Векторне поле $\mathbf{F}$ називається **потенціальним (conservative)** в однозв'язній області, якщо існує така скалярна функція $\Phi(x, y, z)$ (**потенціал**), що:
$$\mathbf{F} = \nabla \Phi \iff P = \frac{\partial \Phi}{\partial x}, \; Q = \frac{\partial \Phi}{\partial y}, \; R = \frac{\partial \Phi}{\partial z}$$

### 3.1. Еквівалентні ознаки потенціальності
1. $\mathbf{F} = \nabla \Phi$.
2. $\operatorname{rot} \mathbf{F} = \mathbf{0}$ (Поле безвихрове).
3. Інтеграл не залежить від траєкторії: $\int_A^B \mathbf{F} \cdot d\mathbf{r} = \Phi(B) - \Phi(A)$ (Аналог формули Ньютона-Лейбніца).
4. Циркуляція по будь-якому замкненому контуру дорівнює нулю: $\oint_C \mathbf{F} \cdot d\mathbf{r} = 0$.

---

## 4. Теорема Гріна (Green's Theorem)

> **👑 Теорема Гріна (на площині $\mathbb{R}^2$):**
> Нехай $D \subset \mathbb{R}^2$ — правильна область, обмежена кусково-гладким контуром $\partial D$, орієнтованим додатно (проти годинникової стрілки). Тоді:
> $$\oint_{\partial D} P\,dx + Q\,dy = \iint_D \left( \frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right) dx\,dy$$

* **Обчислення площі за допомогою інтеграла по межі:**
  $$S(D) = \frac{1}{2}\oint_{\partial D} (x\,dy - y\,dx) = \oint_{\partial D} x\,dy = -\oint_{\partial D} y\,dx$$

---

## 5. Поверхневі інтеграли та потік (Surface Integrals & Flux)

### 5.1. Поверхневий інтеграл 1-го роду (за площею поверхні $dS$)
Для поверхні $S$, заданої рівнянням $z = g(x, y)$ над областю $D$:
$$\iint_S f(x, y, z)\,dS = \iint_D f(x, y, g(x, y)) \sqrt{1 + g_x^2 + g_y^2}\,dx\,dy$$

### 5.2. Поверхневий інтеграл 2-го роду (Потік векторного поля / Flux)
**Потік (Flux $\Phi_{\mathbf{F}}$)** векторного поля $\mathbf{F}$ через орієнтовану поверхню $S$ з одиничним вектором нормалі $\mathbf{n}$:
$$\Phi_{\mathbf{F}} = \iint_S \mathbf{F} \cdot d\mathbf{S} = \iint_S (\mathbf{F} \cdot \mathbf{n})\,dS = \iint_S P\,dy\,dz + Q\,dz\,dx + R\,dx\,dy$$

```text
               z ^       ^ n (Одинична нормаль)
                 |      /
                 |     *-----> F (Векторне поле)
                 |    / \
                 |   / S \ (Орієнтована поверхня)
                 +--'-----'-----------> y
```

---

## 6. Теореми Остроградського-Гаусса та Стокса (Gauss & Stokes)

### 6.1. Теорема Остроградського-Гаусса (Divergence Theorem)
> **👑 Теорема Гаусса-Остроградського (Зв'язок 2D поверхні та 3D об'єму):**
> Потік векторного поля через замкнену поверхню $\partial V$, що обмежує тіло $V$, дорівнює потрійному інтегралу від дивергенції поля по всьому об'єму $V$:
> $$\oiint_{\partial V} \mathbf{F} \cdot d\mathbf{S} = \iiint_V (\operatorname{div} \mathbf{F})\,dV = \iiint_V \left( \frac{\partial P}{\partial x} + \frac{\partial Q}{\partial y} + \frac{\partial R}{\partial z} \right) dV$$

* *Фізичний зміст:* Скільки "рідини" витекло крізь замкнену оболонку, стільки її сумарно нагенерували всі внутрішні джерела.

### 6.2. Теорема Стокса (Stokes' Theorem)
> **👑 Теорема Стокса (Зв'язок 1D контуру та 2D натягнутої поверхні):**
> Циркуляція векторного поля по замкненому просторовому контуру $\partial S$ дорівнює потоку ротора цього поля через будь-яку орієнтовану поверхню $S$, натягнуту на цей контур:
> $$\oint_{\partial S} \mathbf{F} \cdot d\mathbf{r} = \iint_S (\operatorname{rot} \mathbf{F}) \cdot d\mathbf{S} = \iint_S (\nabla \times \mathbf{F}) \cdot \mathbf{n}\,dS$$

* *Зауваження:* Теорема Гріна є окремим 2D випадком теореми Стокса (коли поверхня лежить у площині $xy$).

### 6.3. Єдина формула інтегрального числення (Узагальнена теорема Стокса для дифформ)
Усі 4 фундаментальні теореми (Ньютона-Лейбніца, Гріна, Стокса, Гаусса-Остроградського) є частинними випадками єдиної теореми Картана для диференціальних форм $\omega$:
$$\int_{\partial \Omega} \omega = \int_\Omega d\omega$$

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Орієнтація контуру та нормалі поверхні (Правило правої руки):**
   При застосуванні теореми Стокса орієнтація межі $\partial S$ та напрямок нормалі $\mathbf{n}$ повинні узгоджуватися за правилом свердлика / правої руки. Неправильний вибір normal vector дасть помилку у знаку (-).
2. **Застосування теореми Гаусса до незамкненої поверхні:**
   Формула $\iint_S \mathbf{F} \cdot d\mathbf{S} = \iiint_V \operatorname{div}\mathbf{F}\,dV$ застосовна **виключно до замкнених** поверхонь ($\oiint$). Якщо поверхня відкрита (наприклад, параболоїд без кришки), потрібно або додати основу-кришку, або рахувати напряму.
3. **Потенціальність поля в неоднозв'язних областях:**
   Поле $\mathbf{F} = \left(-\frac{y}{x^2+y^2}, \frac{x}{x^2+y^2}\right)$ має $\operatorname{rot}\mathbf{F} = 0$ всюди крім $(0,0)$, проте інтеграл по одиничному колу $\oint \mathbf{F} \cdot d\mathbf{r} = 2\pi \neq 0$ через наявність особливої точки (дірки) в області.

---

[⬅️ Попередня: 08. Кратні інтеграли та заміни змінних](08-multiple-integrals-and-coordinate-transformations.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Диференціальні рівняння та системи ➡️](10-differential-equations-ode-and-systems.md) | [⚡ Cheat Sheet](cheat-sheet.md)
