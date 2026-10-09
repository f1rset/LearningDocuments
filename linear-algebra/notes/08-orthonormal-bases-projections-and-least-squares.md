# 08. Ортонормовані базиси, ортогональні проекції та метод найменших квадратів (Orthonormal Bases, Projections & Least Squares)

[⬅️ 07. Норми та ортогональність](07-norms-inner-products-and-orthogonality.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 09. Ортогональні матриці, Грама-Шмідта та QR-розклад ➡️](09-orthogonal-matrices-gram-schmidt-and-qr.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Ортонормований базис (Orthonormal Basis / ONB):** Базис із взаємно перпендикулярних одиничних векторів.
2. **Ортогональна проекція (Orthogonal Projection):** Найкраще наближення вектора у підпросторі.
3. **Матриця проекції (Projection Matrix $P$):** Властивості $P^2 = P$ (ідемпотентність) та $P^T = P$ (симетричність).
4. **Метод найменших квадратів (Least Squares Method / OLS):** Знаходження найкращого псевдорозв'язку для перевизначених систем $A\mathbf{x} \approx \mathbf{b}$.
5. **Нормальні рівняння (Normal Equations):** $A^T A \hat{\mathbf{x}} = A^T \mathbf{b}$.

---

## 1. Ортонормовані базиси (Orthonormal Bases)

### 1.1. Означення
Множина векторів $\{\mathbf{q}_1, \mathbf{q}_2, \dots, \mathbf{q}_n\} \subset V$ називається **ортонормованою (orthonormal)**, якщо всі вектори попарно ортогональні та мають одиничну довжину (нормовані):
$$\langle \mathbf{q}_i, \mathbf{q}_j \rangle = \mathbf{q}_i^T \mathbf{q}_j = \delta_{ij} = \begin{cases} 1, & \text{якщо } i = j \\ 0, & \text{якщо } i \neq j \end{cases}$$
де $\delta_{ij}$ — символ Кронекера (Kronecker delta).

### 1.2. Чому ортонормовані базиси є найкращими? (Суперсила ONB)
Для довільного базису знаходження координат $[\mathbf{x}]_{\mathcal{B}}$ вимагає розв'язання системи рівнянь ($P^{-1} \mathbf{x}$, складність $\mathcal{O}(n^3)$).
Для **ортонормованого базису** координати розраховуються **миттєво через скалярні добутки**:
$$\mathbf{x} = \sum_{i=1}^n c_i \mathbf{q}_i \implies \mathbf{c_i = \langle \mathbf{x}, \mathbf{q}_i \rangle = \mathbf{q}_i^T \mathbf{x}}$$
$$\mathbf{x} = (\mathbf{q}_1^T \mathbf{x})\mathbf{q}_1 + (\mathbf{q}_2^T \mathbf{x})\mathbf{q}_2 + \dots + (\mathbf{q}_n^T \mathbf{x})\mathbf{q}_n$$

* **Рівність Парсеваля (Parseval's Identity):**
  $$\|\mathbf{x}\|^2 = \sum_{i=1}^n c_i^2 = c_1^2 + c_2^2 + \dots + c_n^2$$

---

## 2. Ортогональні проекції (Orthogonal Projections)

### 2.1. Проекція вектора на пряму (одновимірний підпростір)
Нехай вектор $\mathbf{a} \neq \mathbf{0}$ задає напрямок прямої. Проекція вектора $\mathbf{b}$ на пряму $\text{span}\{\mathbf{a}\}$:

```text
               b ^
                 |  \
                 |    \  Помилка e = b - p (e _|_ a)
                 |      \
                 +-------+---------> a
                 0       p = x_hat * a
```

1. Вектор проекції $\mathbf{p}$ лежить на прямій: $\mathbf{p} = \hat{x} \mathbf{a}$.
2. Вектор похибки (помилки) $\mathbf{e} = \mathbf{b} - \mathbf{p}$ ортогональний до $\mathbf{a}$:
   $$\mathbf{a}^T (\mathbf{b} - \hat{x} \mathbf{a}) = 0 \iff \mathbf{a}^T \mathbf{b} - \hat{x} \mathbf{a}^T \mathbf{a} = 0 \implies \hat{x} = \frac{\mathbf{a}^T \mathbf{b}}{\mathbf{a}^T \mathbf{a}} = \frac{\mathbf{a}^T \mathbf{b}}{\|\mathbf{a}\|^2}$$
3. Вектор проекції:
   $$\mathbf{p} = \hat{x} \mathbf{a} = \frac{\mathbf{a}^T \mathbf{b}}{\mathbf{a}^T \mathbf{a}} \mathbf{a} = \left( \frac{\mathbf{a} \mathbf{a}^T}{\mathbf{a}^T \mathbf{a}} \right) \mathbf{b} = P \mathbf{b}$$
4. **Матриця проекції на пряму:**
   $$P = \frac{\mathbf{a} \mathbf{a}^T}{\mathbf{a}^T \mathbf{a}}$$

---

### 2.2. Проекція вектора на довільний підпростір $W = \text{Col}(A)$
Нехай підпростір $W \le \mathbb{R}^m$ натягнутий на лінійно незалежні стовпці матриці $A \in \mathbb{R}^{m \times n}$ ($m > n$).
Ми шукаємо проекцію $\mathbf{p} \in \text{Col}(A)$, найближчу до $\mathbf{b} \in \mathbb{R}^m$.

1. Оскільки $\mathbf{p} \in \text{Col}(A)$, існує $\hat{\mathbf{x}} \in \mathbb{R}^n$ такий, що $\mathbf{p} = A \hat{\mathbf{x}}$.
2. Вектор похибки $\mathbf{e} = \mathbf{b} - A\hat{\mathbf{x}}$ має бути **перпендикулярним до всього підпростору $\text{Col}(A)$**, тобто до кожного стовпця матриці $A$:
   $$A^T \mathbf{e} = \mathbf{0} \iff A^T (\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}$$
3. Звідси отримуємо **Нормальні рівняння (Normal Equations)**:
   $$\mathbf{A^T A \hat{\mathbf{x}} = A^T \mathbf{b}}$$
4. Оскільки стовпці $A$ лінійно незалежні, симетрична квадратна матриця $A^T A \in \mathbb{R}^{n \times n}$ є оборотною ($(A^TA)^{-1}$ існує):
   $$\mathbf{\hat{\mathbf{x}} = (A^T A)^{-1} A^T \mathbf{b}}$$
5. Вектор проекції:
   $$\mathbf{p} = A \hat{\mathbf{x}} = \mathbf{A (A^T A)^{-1} A^T \mathbf{b}} = P \mathbf{b}$$
6. **Матриця проекції (Projection Matrix $P$ на підпростір $\text{Col}(A)$):**
   $$\mathbf{P = A (A^T A)^{-1} A^T}$$

### 2.3. Властивості проекційної матриці $P$:
1. **Симетричність (Symmetric):** $P^T = (A(A^TA)^{-1}A^T)^T = A((A^TA)^{-1})^T A^T = P$.
2. **Ідемпотентність (Idempotent):** Повторна проекція не змінює результат:
   $$P^2 = P \cdot P = [A(A^TA)^{-1}A^T] [A(A^TA)^{-1}A^T] = A(A^TA)^{-1} (A^TA) (A^TA)^{-1} A^T = P$$
3. **Проекція на ортогональне доповнення $W^\perp$:**
   $$P_\perp = I - P \implies \mathbf{e} = (I - P)\mathbf{b}$$
4. Власні значення проекційної матриці завжди дорівнюють або **$0$**, або **$1$**: $\lambda \in \{0, 1\}$.

---

## 3. Метод найменших квадратів (Least Squares Approximation / OLS)

### 3.1. Постановка задачі
Нехай система $A\mathbf{x} = \mathbf{b}$ є **перевизначеною (overdetermined)** ($m > n$, кількість експериментальних спостережень $m$ більша за кількість невідомих параметрів $n$) і **несумісною** ($\mathbf{b} \notin \text{Col}(A)$). Точного розв'язку не існує.

**Задача найменших квадратів:** Знайти вектор $\hat{\mathbf{x}}$, який **мінімізує норму нев'язки (квадрат похибки)**:
$$\min_{\mathbf{x} \in \mathbb{R}^n} \|A\mathbf{x} - \mathbf{b}\|_2^2 = \min_{\mathbf{x}} \sum_{i=1}^m ( (A\mathbf{x})_i - b_i )^2$$

```text
       Геометрія методу найменших квадратів:
       
             b (Вектор реальних вимірювань R^m)
             ^
             | \
             |   \  Похибка e = b - p (МІНІМАЛЬНА, e _|_ Col(A))
             |     \
         ----+-------+--------------------> Col(A) (Підпростір моделей)
            (0)      p = A * x_hat (Найкраща проекція)
```

### 3.2. Теорема про найкраще наближення (Best Approximation Theorem)
Вектор $\hat{\mathbf{x}}$ є розв'язком задачі найменших квадратів тоді і тільки тоді, коли він є розв'язком системи **нормальних рівнянь**:
$$A^T A \hat{\mathbf{x}} = A^T \mathbf{b} \iff \hat{\mathbf{x}} = (A^T A)^{-1} A^T \mathbf{b} = A^+ \mathbf{b}$$
де $A^+ = (A^T A)^{-1} A^T$ — ліва **псевдообернена матриця Мура-Пенроуза (Moore-Penrose Pseudoinverse)**.

### 3.3. Приклад: Лінійна регресія (Linear Regression / Line Fitting)
Підібрати пряму $y = \beta_0 + \beta_1 x$ для набору точок $(x_1, y_1), \dots, (x_m, y_m)$:
$$
\begin{bmatrix} 1 & x_1 \\ 1 & x_2 \\ \vdots & \vdots \\ 1 & x_m \end{bmatrix} \begin{bmatrix} \beta_0 \\ \beta_1 \end{bmatrix} \approx \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_m \end{bmatrix} \iff A \boldsymbol{\beta} \approx \mathbf{y} \implies \boldsymbol{\beta} = (A^T A)^{-1} A^T \mathbf{y}
$$

---

## 4. ⚠️ Типові підводні камені (Common Pitfalls)
* 🚨 **Спроба скоротити $(A^T A)^{-1}$ як $A^{-1} (A^T)^{-1}$:** Ця формула **НЕПРАВИЛЬНА** для прямокутних матриць ($m \neq n$), оскільки для прямокутної матриці $A$ окрема обернена $A^{-1}$ взагалі не існує! Матриця $A^T A$ квадратна ($n \times n$), і обертати можна тільки весь добуток разом.

---

## 📊 Зведена таблиця (Summary Table)

| Поняття | Формула | Властивість |
|---|---|---|
| Ортонормований базис | $\mathbf{q}_i^T \mathbf{q}_j = \delta_{ij}$ | $\mathbf{x} = \sum (\mathbf{q}_i^T \mathbf{x}) \mathbf{q}_i$ |
| Проекція на пряму | $\mathbf{p} = \frac{\mathbf{a}^T \mathbf{b}}{\mathbf{a}^T \mathbf{a}} \mathbf{a}$ | $P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T \mathbf{a}}$ |
| Нормальні рівняння | $A^T A \hat{\mathbf{x}} = A^T \mathbf{b}$ | Мінімізує $\|A\mathbf{x} - \mathbf{b}\|_2^2$ |
| Матриця проекції на $\text{Col}(A)$ | $P = A(A^TA)^{-1}A^T$ | $P^2 = P, \quad P^T = P$ |
| Псевдообернена Мура-Пенроуза | $A^+ = (A^T A)^{-1} A^T$ | $\hat{\mathbf{x}} = A^+ \mathbf{b}$ |

---

[⬅️ 07. Норми та ортогональність](07-norms-inner-products-and-orthogonality.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 09. Ортогональні матриці, Грама-Шмідта та QR-розклад ➡️](09-orthogonal-matrices-gram-schmidt-and-qr.md) | [⚡ Cheat Sheet](cheat-sheet.md)
