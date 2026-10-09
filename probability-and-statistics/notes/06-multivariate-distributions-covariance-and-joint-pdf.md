# 06. Багатовимірні розподіли, коваріація та спільна щільність (Multivariate Distributions, Covariance & Joint PDF)

[⬅️ Попередня: 05. Числові характеристики та моменти](05-expectation-variance-and-moments.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Граничні теореми: Закон великих чисел та ЦГТ ➡️](07-limit-theorems-law-of-large-numbers-and-clt.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Спільні розподіли (Joint PMF / Joint PDF $f_{X,Y}(x, y)$):** Імовірнісний опис взаємопов'язаних випадкових величин.
2. **Маргінальні розподіли (Marginal Distributions):** Проєктування спільного розподілу на окремі змінні.
3. **Умовний розподіл та закон повного сподівання (Law of Total Expectation):** $E[E[Y|X]] = E[Y]$.
4. **Коваріація та коефіцієнт кореляції Пірсона ($\rho \in [-1, 1]$):** Лінійна залежність та відмінність від незалежності.
5. **Коваріаційна матриця ($\Sigma \succeq 0$):** Багатовимірне розсіювання та дисперсія лінійних комбінацій $\operatorname{Var}(\mathbf{a}^T\mathbf{X}) = \mathbf{a}^T\Sigma\mathbf{a}$.
6. **Багатовимірний нормальний розподіл ($\mathcal{N}(\boldsymbol{\mu}, \Sigma)$):** Геометрія еліпсоїдів розсіювання та еквівалентність некорельованості й незалежності.

---

## 1. Спільні та маргінальні розподіли (Joint & Marginal Distributions)

### 1.1. Спільна щільність (Joint PDF)
Для випадкового вектора $(X, Y) \in \mathbb{R}^2$ існує функція $f_{X,Y}(x, y)$, така що:
$$P((X, Y) \in A) = \iint_A f_{X,Y}(x, y)\,dx\,dy$$
* **Умова нормування:** $\int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} f_{X,Y}(x, y)\,dx\,dy = 1$.
* **Спільна функція розподілу (Joint CDF):**
  $$F_{X,Y}(x, y) = P(X \le x, Y \le y) = \int_{-\infty}^x \int_{-\infty}^y f_{X,Y}(u, v)\,du\,dv$$

### 1.2. Маргінальні розподіли (Marginal Distributions)
Отримання розподілу однієї змінної шляхом інтегрування (або сумування) по іншій:
* **Неперервний випадок:**
  $$f_X(x) = \int_{-\infty}^{+\infty} f_{X,Y}(x, y)\,dy, \qquad f_Y(y) = \int_{-\infty}^{+\infty} f_{X,Y}(x, y)\,dx$$
* **Дискретний випадок:**
  $$p_X(x) = \sum_y p_{X,Y}(x, y), \qquad p_Y(y) = \sum_x p_{X,Y}(x, y)$$

---

## 2. Умовні розподіли та умовне математичне сподівання

### 2.1. Умовна щільність (Conditional PDF)
$$f_{Y \mid X}(y \mid x) = \frac{f_{X,Y}(x, y)}{f_X(x)}, \quad \text{за умови } f_X(x) > 0$$

### 2.2. Умовне математичне сподівання (Conditional Expectation)
$E[Y \mid X = x] = \int_{-\infty}^\infty y f_{Y|X}(y \mid x)\,dy$ є числом, а $E[Y \mid X]$ є **випадковою величиною** (функцією від $X$).

* **👑 Закон повного математичного сподівання (Law of Total Expectation / Adam's Law):**
  $$E\left[ E[Y \mid X] \right] = E[Y]$$
* **👑 Закон повної дисперсії (Law of Total Variance / Eve's Law):**
  $$\operatorname{Var}(Y) = E\left[ \operatorname{Var}(Y \mid X) \right] + \operatorname{Var}\left( E[Y \mid X] \right)$$
  *(Дисперсія = Середнє внутрішньогрупової дисперсії + Дисперсія міжгрупових середніх)*.

---

## 3. Незалежність випадкових величин (Independence)

Випадкові величини $X$ та $Y$ називаються **незалежними**, якщо їх спільна щільність факторизується в добуток маргінальних:
$$f_{X,Y}(x, y) = f_X(x) \cdot f_Y(y) \iff F_{X,Y}(x, y) = F_X(x) \cdot F_Y(y), \quad \forall (x, y)$$

---

## 4. Коваріація та кореляція (Covariance & Correlation)

### 4.1. Коваріація (Covariance)
Міра лінійної спільної мінливості двох величин:
$$\operatorname{Cov}(X, Y) = \sigma_{XY} = E\left[ (X - \mu_X)(Y - \mu_Y) \right] = E[X Y] - E[X] E[Y]$$

* **Властивості коваріації:**
  1. $\operatorname{Cov}(X, X) = \operatorname{Var}(X)$.
  2. Симетрія: $\operatorname{Cov}(X, Y) = \operatorname{Cov}(Y, X)$.
  3. Білінійність: $\operatorname{Cov}(aX + b, cY + d) = a c \operatorname{Cov}(X, Y)$.
  4. Якщо $X$ та $Y$ незалежні $\implies \operatorname{Cov}(X, Y) = 0$.

### 4.2. Коефіцієнт лінійної кореляції Пірсона (Pearson Correlation Coefficient $\rho$)
$$\rho(X, Y) = \operatorname{Corr}(X, Y) = \frac{\operatorname{Cov}(X, Y)}{\sigma_X \sigma_Y}$$

* **Властивості кореляції:**
  1. $-1 \le \rho(X, Y) \le 1$ (випливає з нерівності Коші-Буняковського).
  2. $|\rho| = 1 \iff Y = a X + b$ ($a > 0 \implies \rho = +1$, $a < 0 \implies \rho = -1$) — строга лінійна залежність.
  3. $\rho = 0 \implies$ величини називаються **некорельованими (uncorrelated)**.

> **⚠️ Критичне розрізнення (Некорельованість $\neq$ Незалежність):**
> * Незалежність $\implies \rho = 0$.
> * Зворотне в загальному випадку **НЕ ВІРНО**! Кореляція вимірює лише **лінійний** зв'язок.
> * *Контрприклад:* Нехай $X \sim U(-1, 1)$, $Y = X^2$. Величини жорстко функціонально залежні, але $\operatorname{Cov}(X, Y) = E[X^3] - E[X]E[X^2] = 0 \implies \rho = 0$.

```text
    ρ = +1 (Пряма лінія)      ρ = 0 (Некорельовані)     ρ = 0 (Нелінійна залежність Y=X^2)
        y ^    /                  y ^  . * . .              y ^  \     /
          |   /                     | * . * . *               |   \   /
          |  /                      |  . * .                  |    \_/
          +-------> x               +-----------> x           +-----------> x
```

---

## 5. Багатовимірний вектор та коваріаційна матриця (Covariance Matrix)

Для випадкового вектора $\mathbf{X} = (X_1, X_2, \dots, X_n)^T \in \mathbb{R}^n$:
* **Вектор середніх:** $\boldsymbol{\mu} = E[\mathbf{X}] = (\mu_1, \dots, \mu_n)^T$.
* **Коваріаційна матриця (Covariance Matrix $\Sigma \in \mathbb{R}^{n \times n}$):**
  $$\Sigma = E\left[ (\mathbf{X} - \boldsymbol{\mu})(\mathbf{X} - \boldsymbol{\mu})^T \right] = \begin{pmatrix} 
  \operatorname{Var}(X_1) & \operatorname{Cov}(X_1, X_2) & \dots & \operatorname{Cov}(X_1, X_n) \\ 
  \operatorname{Cov}(X_2, X_1) & \operatorname{Var}(X_2) & \dots & \operatorname{Cov}(X_2, X_n) \\ 
  \vdots & \vdots & \ddots & \vdots \\ 
  \operatorname{Cov}(X_n, X_1) & \operatorname{Cov}(X_n, X_2) & \dots & \operatorname{Var}(X_n) 
  \end{pmatrix}$$

* **Властивості матриці $\Sigma$:**
  1. Симетрична: $\Sigma = \Sigma^T$.
  2. Додатно напіввизначена: $\Sigma \succeq 0$ (усі власні значення $\lambda_i \ge 0$).
  3. **Дисперсія лінійної комбінації:** Для будь-якого вектора ваг $\mathbf{a} \in \mathbb{R}^n$:
     $$\operatorname{Var}(\mathbf{a}^T \mathbf{X}) = \mathbf{a}^T \Sigma \mathbf{a} \ge 0$$

---

## 6. Багатовимірний нормальний розподіл (Multivariate Gaussian $\mathcal{N}(\boldsymbol{\mu}, \Sigma)$)

Вектор $\mathbf{X} \sim \mathcal{N}_n(\boldsymbol{\mu}, \Sigma)$ має спільну щільність:
$$f(\mathbf{x}) = \frac{1}{(2\pi)^{n/2} |\det\Sigma|^{1/2}} \exp\left( -\frac{1}{2} (\mathbf{x} - \boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x} - \boldsymbol{\mu}) \right)$$

* **Властивості Гауссового вектора:**
  1. **Лінійне перетворення:** Якщо $\mathbf{Y} = A \mathbf{X} + \mathbf{b}$, то $\mathbf{Y} \sim \mathcal{N}(A\boldsymbol{\mu} + \mathbf{b}, \; A \Sigma A^T)$.
  2. **👑 Виняток для Гауссових величин:** Якщо вектор $\mathbf{X}$ є сумісно нормальним, то **некорельованість компонент еквівалентна їх повній незалежності** ($\Sigma$ діагональна $\iff X_i$ незалежні)!
  3. **Еліпсоїди Махаланобіса:** Лінії однакової щільності утворюють еліпсоїди $(\mathbf{x} - \boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x} - \boldsymbol{\mu}) = c^2$, головні осі яких збігаються з власними векторами матриці $\Sigma$.

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Маргінальна нормальність НЕ гарантує спільну нормальність:**
   Якщо кожна величина $X$ та $Y$ окремо має нормальний розподіл, це ще **не означає**, що їх спільний вектор $(X, Y)$ є 2D-нормальним.
2. **Невід'ємність коваріаційної матриці:** Довільна таблиця чисел не може бути коваріаційною матрицею — вона зобов'язана бути симетричною і додатно напіввизначеною ($\det \Sigma \ge 0$).
3. **Плутанина між формулами дисперсії суми:**
   $\operatorname{Var}(X - Y) = \operatorname{Var}(X) + \operatorname{Var}(Y) - 2\operatorname{Cov}(X, Y)$ (знак мінус перед коваріацією, але плюс між дисперсіями!).

---

[⬅️ Попередня: 05. Числові характеристики та моменти](05-expectation-variance-and-moments.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Граничні теореми: Закон великих чисел та ЦГТ ➡️](07-limit-theorems-law-of-large-numbers-and-clt.md) | [⚡ Cheat Sheet](cheat-sheet.md)

