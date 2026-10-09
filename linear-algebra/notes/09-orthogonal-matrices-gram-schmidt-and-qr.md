# 09. Ортогональні матриці, процес Грама-Шмідта та QR-факторизація (Orthogonal Matrices, Gram-Schmidt & QR Factorization)

[⬅️ 08. Проекції та метод найменших квадратів](08-orthonormal-bases-projections-and-least-squares.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Власні значення, вектори та жорданова форма ➡️](10-eigenvalues-eigenvectors-and-jordan-form.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Ортогональна матриця (Orthogonal Matrix $Q$):** Квадратна матриця з ортонормованими стовпцями ($Q^T Q = Q Q^T = I$).
2. **Збереження довжин та кутів (Isometry):** $\|Q\mathbf{x}\| = \|\mathbf{x}\|$ та $\langle Q\mathbf{x}, Q\mathbf{y} \rangle = \langle \mathbf{x}, \mathbf{y} \rangle$.
3. **Процес ортогоналізації Грама-Шмідта (Gram-Schmidt Orthogonalization Process):** Алгоритм перетворення довільного базису в ортонормований.
4. **QR-розклад (QR Factorization):** Розклад $A = QR$ (ортогональна матриця $\times$ верхня трикутна).
5. **Чисельно стійкий МНК (Stable Least Squares via QR):** $R\hat{\mathbf{x}} = Q^T \mathbf{b}$.

---

## 1. Ортогональні матриці (Orthogonal Matrices)

### 1.1. Означення
Квадратна дійсна матриця $Q \in \mathbb{R}^{n \times n}$ називається **ортогональною (orthogonal matrix)**, якщо її стовпці утворюють ортонормований базис в $\mathbb{R}^n$:
$$\mathbf{Q^T Q = I_n \iff Q^{-1} = Q^T}$$
*(Для комплексних просторів аналогом є **унітарна матриця (Unitary Matrix $U$)**: $U^* U = I$, де $U^* = \overline{U}^T$).*

### 1.2. Фундаментальні властивості ортогональних матриць:
1. **Легкість обернення:** Обернена матриця дорівнює транспонованій: $Q^{-1} = Q^T$ (знаходиться за $\mathcal{O}(1)$ операцій без жодного ділення!).
2. **Стовпці та рядки:** І стовпці, і рядки матриці $Q$ є ортонормованими.
3. **Збереження скалярного добутку:**
   $$\langle Q\mathbf{x}, Q\mathbf{y} \rangle = (Q\mathbf{x})^T (Q\mathbf{y}) = \mathbf{x}^T Q^T Q \mathbf{y} = \mathbf{x}^T I \mathbf{y} = \mathbf{x}^T \mathbf{y} = \langle \mathbf{x}, \mathbf{y} \rangle$$
4. **Збереження довжини (Ізометрія / Isometry):**
   $$\|Q\mathbf{x}\|_2 = \sqrt{(Q\mathbf{x})^T (Q\mathbf{x})} = \sqrt{\mathbf{x}^T \mathbf{x}} = \|\mathbf{x}\|_2$$
5. **Збереження кутів та відстаней:** $d(Q\mathbf{x}, Q\mathbf{y}) = \|\mathbf{x} - \mathbf{y}\|_2$.
6. **Визначник:** $\det(Q^T Q) = \det(I) \implies (\det(Q))^2 = 1 \implies \mathbf{\det(Q) = \pm 1}$.
   * $\det(Q) = +1$: **Чисте обертання (Rotation)** (спеціальна ортогональна група $SO(n)$).
   * $\det(Q) = -1$: **Відбиття (Reflection / Mirroring)** (або обертання з відбиттям).

---

## 2. Процес ортогоналізації Грама-Шмідта (Gram-Schmidt Process)

### 2.1. Мета алгоритму
Перетворити довільний набір лінійно незалежних векторів $\{\mathbf{a}_1, \mathbf{a}_2, \dots, \mathbf{a}_n\}$ у **взаємно перпендикулярний ортонормований набір** $\{\mathbf{q}_1, \mathbf{q}_2, \dots, \mathbf{q}_n\}$ із збереженням підпросторів: $\text{span}\{\mathbf{q}_1, \dots, \mathbf{q}_k\} = \text{span}\{\mathbf{a}_1, \dots, \mathbf{a}_k\}$.

```text
       Геометрія кроку Грама-Шмідта для 2 векторів:
       
             a2 ^
                | \
                |   \  u2 = a2 - p1 (Ортогональний залишок, u2 _|_ q1)
                |     \
         -------+-------+-----------------> a1 || q1
               (0)      p1 = (q1^T * a2) * q1 (Проекція a2 на q1)
```

### 2.2. Покроковий алгоритм (Classical Gram-Schmidt / CGS):

* **Крок 1:** Перший напрямок беремо без змін та нормуємо:
  $$\mathbf{u}_1 = \mathbf{a}_1, \quad \mathbf{q}_1 = \frac{\mathbf{u}_1}{\|\mathbf{u}_1\|}$$

* **Крок 2:** Від вектора $\mathbf{a}_2$ віднімаємо його проекцію на $\mathbf{q}_1$:
  $$\mathbf{u}_2 = \mathbf{a}_2 - (\mathbf{q}_1^T \mathbf{a}_2) \mathbf{q}_1, \quad \mathbf{q}_2 = \frac{\mathbf{u}_2}{\|\mathbf{u}_2\|}$$

* **Крок 3:** Від вектора $\mathbf{a}_3$ віднімаємо його проекції на вже знайдені $\mathbf{q}_1$ та $\mathbf{q}_2$:
  $$\mathbf{u}_3 = \mathbf{a}_3 - (\mathbf{q}_1^T \mathbf{a}_3) \mathbf{q}_1 - (\mathbf{q}_2^T \mathbf{a}_3) \mathbf{q}_2, \quad \mathbf{q}_3 = \frac{\mathbf{u}_3}{\|\mathbf{u}_3\|}$$

* **Крок $k$ (Загальна формула):**
  $$\mathbf{u}_k = \mathbf{a}_k - \sum_{i=1}^{k-1} (\mathbf{q}_i^T \mathbf{a}_k) \mathbf{q}_i, \quad \mathbf{q}_k = \frac{\mathbf{u}_k}{\|\mathbf{u}_k\|}$$

> **💡 Модифікований алгоритм Грама-Шмідта (Modified Gram-Schmidt / MGS):**
> У комп'ютерній floating-point арифметиці класичний алгоритм CGS швидко втрачає ортогональність через накопичення похибок заокруглення. Алгоритм MGS проектує вектори послідовно на кожен щойно оновлений залишок, що забезпечує **високу чисельну стійкість**.

---

## 3. QR-розклад (QR Factorization)

### 3.1. Матрична форма Грама-Шмідта
Виразимо вихідні вектори $\mathbf{a}_k$ через знайдені ортонормовані вектори $\mathbf{q}_i$:
$$\mathbf{a}_1 = (\mathbf{q}_1^T \mathbf{a}_1) \mathbf{q}_1 = r_{11} \mathbf{q}_1$$
$$\mathbf{a}_2 = (\mathbf{q}_1^T \mathbf{a}_2) \mathbf{q}_1 + (\mathbf{q}_2^T \mathbf{a}_2) \mathbf{q}_2 = r_{12} \mathbf{q}_1 + r_{22} \mathbf{q}_2$$
$$\mathbf{a}_k = \sum_{i=1}^k r_{ik} \mathbf{q}_i, \quad \text{де } r_{ik} = \mathbf{q}_i^T \mathbf{a}_k \ (i < k), \quad r_{kk} = \|\mathbf{u}_k\|$$

У матричному вигляді для $A \in \mathbb{R}^{m \times n}$ ($m \ge n$):
$$\mathbf{A = Q R}$$
де:
* $Q \in \mathbb{R}^{m \times n}$ — матриця з **ортонормованими стовпцями** ($Q^T Q = I_n$).
* $R \in \mathbb{R}^{n \times n}$ — **верхня трикутна матриця (Upper triangular matrix)** з додатними елементами на діагоналі ($r_{kk} > 0$):
  $$R = \begin{bmatrix} r_{11} & r_{12} & \dots & r_{1n} \\ 0 & r_{22} & \dots & r_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & r_{nn} \end{bmatrix}$$

---

## 4. Застосування QR-розкладу для розв'язання МНК (Least Squares via QR)

### 4.1. Чому МНК через QR кращий за нормальні рівняння?
Звичайні нормальні рівняння вимагають обчислення матриці $A^T A$:
* **Проблема обумовленості (Condition Number):** Число обумовленості матриці підноситься до квадрата: $\kappa(A^T A) = (\kappa(A))^2$. Якщо матриця $A$ була погано обумовленою, $A^T A$ втрачає половину значущих цифр точності!
* **Розв'язання через QR:**
  Підставимо $A = QR$ у нормальні рівняння:
  $$A^T A \hat{\mathbf{x}} = A^T \mathbf{b} \iff (QR)^T (QR) \hat{\mathbf{x}} = (QR)^T \mathbf{b}$$
  $$R^T \underbrace{Q^T Q}_{I} R \hat{\mathbf{x}} = R^T Q^T \mathbf{b} \iff R^T R \hat{\mathbf{x}} = R^T Q^T \mathbf{b}$$
  Оскільки $R$ невироджена, домножимо зліва на $(R^T)^{-1}$:
  $$\mathbf{R \hat{\mathbf{x}} = Q^T \mathbf{b}}$$

> **🔥 ПЕРЕВАГА:**
> Замість розв'язання складної системи з $A^T A$, ми просто:
> 1. Множимо $\mathbf{d} = Q^T \mathbf{b}$ (ортогональна проекція).
> 2. Розв'язуємо трикутну систему $R \hat{\mathbf{x}} = \mathbf{d}$ **миттєвою зворотною підстановкою (back substitution)**!

---

## 5. ⚠️ Типові підводні камені
* 🚨 **Прямокутна матриця $Q$ ($m > n$):** Якщо $Q \in \mathbb{R}^{m \times n}$, то $Q^T Q = I_n$ (одинична $n \times n$), але $Q Q^T \neq I_m$! Матриця $P = Q Q^T \in \mathbb{R}^{m \times m}$ є **матрицею проекції** на підпростір $\text{Col}(A)$.

---

## 📊 Зведена таблиця (Summary Table)

| Властивість / Розклад | Формула | Схемотехнічне / Алгоритмічне значення |
|---|---|---|
| Ортогональна матриця | $Q^T Q = I \iff Q^{-1} = Q^T$ | Зберігає довжини, обернення без обчислень |
| Грама-Шмідта | $\mathbf{u}_k = \mathbf{a}_k - \sum (\mathbf{q}_i^T \mathbf{a}_k)\mathbf{q}_i$ | Створення ортонормованого базису |
| QR-розклад | $A = QR$ | $Q^T Q = I, \quad R$ — upper triangular |
| МНК через QR | $R\hat{\mathbf{x}} = Q^T \mathbf{b}$ | Найбільш чисельно стійкий алгоритм лінійної регресії |
| Матриця проекції через $Q$ | $P = Q Q^T$ | Спрощення $A(A^TA)^{-1}A^T \to Q Q^T$ |

---

[⬅️ 08. Проекції та метод найменших квадратів](08-orthonormal-bases-projections-and-least-squares.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Власні значення, вектори та жорданова форма ➡️](10-eigenvalues-eigenvectors-and-jordan-form.md) | [⚡ Cheat Sheet](cheat-sheet.md)

