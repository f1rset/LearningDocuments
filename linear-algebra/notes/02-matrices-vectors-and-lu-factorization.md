# 02. Матриці, вектори, елементарні перетворення та LU-розклад (Matrices, Vectors, Elementary Transformations & LU Factorization)

[⬅️ 01. Системи лінійних рівнянь](01-systems-of-linear-equations-and-gaussian-elimination.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 03. Визначники ➡️](03-determinants-and-cramers-rule.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Матричні операції (Matrix Operations):** Додавання, множення на скаляр, матричне множення, транспонування.
2. **Елементарні матриці (Elementary Matrices):** Матричне представлення елементарних рядкових перетворень.
3. **Оборотні та необоротні матриці (Invertible & Non-invertible / Singular Matrices):** Критерії оборотності.
4. **LU-розклад (LU Factorization) та LUP-розклад (LUP with Row Permutation):** Розклад матриці на нижню ($L$) та верхню ($U$) трикутні матриці.

---

## 1. Матриці та вектори: базові операції (Matrices & Vectors)

### 1.1. Означення та розмірності
* Матриця $A \in \mathbb{R}^{m \times n}$ складається з $m$ рядків та $n$ стовпців:
  $$A = [a_{ij}] = \begin{bmatrix} a_{11} & a_{12} & \dots & a_{1n} \\ a_{21} & a_{22} & \dots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \dots & a_{mn} \end{bmatrix}$$
* Вектор-стовпець (column vector) $\mathbf{x} \in \mathbb{R}^n \equiv \mathbb{R}^{n \times 1}$, вектор-рядок (row vector) $\mathbf{x}^T \in \mathbb{R}^{1 \times n}$.

### 1.2. Основні операції над матрицями
1. **Лінійні операції:** $A + B = [a_{ij} + b_{ij}]$, $\alpha A = [\alpha a_{ij}]$.
2. **Добуток матриць (Matrix Multiplication):**
   Нехай $A \in \mathbb{R}^{m \times k}$ та $B \in \mathbb{R}^{k \times n}$. Тоді $C = AB \in \mathbb{R}^{m \times n}$:
   $$c_{ij} = \sum_{p=1}^k a_{ip} b_{pj} = \mathbf{a}_{i,:} \cdot \mathbf{b}_{:,j}$$
   * **4 погляди на множення матриць:**
     1. **Елементний (Entry):** $c_{ij}$ — скалярний добуток $i$-го рядка $A$ на $j$-й стовпець $B$.
     2. **Стовпчиковий (Columns of $AB$):** Кожен стовпець $C_{:,j} = A \mathbf{b}_{:,j}$ — лінійна комбінація стовпців $A$.
     3. **Рядковий (Rows of $AB$):** Кожен рядок $C_{i,:} = \mathbf{a}_{i,:} B$ — лінійна комбінація рядків $B$.
     4. **Зовнішній добуток (Outer Product):** $AB = \sum_{p=1}^k \mathbf{a}_{:,p} \mathbf{b}_{p,:}$ — сума матриць рангу 1.
   * **Властивості:** Асоціативність $A(BC) = (AB)C$, дистрибутивність $A(B+C) = AB + AC$, але **некомутативність у загальному випадку:** $AB \neq BA$.
3. **Транспонування (Transpose):** $A^T \in \mathbb{R}^{n \times m}$, де $(A^T)_{ij} = a_{ji}$.
   * Властивості: $(A + B)^T = A^T + B^T$, $(AB)^T = B^T A^T$, $(A^T)^T = A$.
   * **Симетрична матриця (Symmetric matrix):** $A^T = A$.
   * **Кососиметрична матриця (Anti-symmetric / Skew-symmetric matrix):** $A^T = -A$.
4. **Слід матриці (Trace):** Сума діагональних елементів квадратної матриці: $\text{tr}(A) = \sum_{i=1}^n a_{ii}$.
   * Властивість циклічності: $\text{tr}(AB) = \text{tr}(BA)$.

---

## 2. Елементарні матриці та зворотна матриця (Elementary & Invertible Matrices)

### 2.1. Елементарні матриці (Elementary Matrices)
**Елементарна матриця $E$** — це квадратна матриця, отримана застосуванням одного елементарного рядкового перетворення до одиничної матриці $I_m$.
* Виконання рядкового перетворення над $A$ еквівалентне множенню **зліва**: $A \to EA$.
* Кожна елементарна матриця є **оборотною (invertible)**, а $E^{-1}$ — елементарна матриця того ж типу:
  1. $E_{i \leftrightarrow j}^{-1} = E_{i \leftrightarrow j}$ (перестановка).
  2. $E_{i \leftarrow c R_i}^{-1} = E_{i \leftarrow \frac{1}{c} R_i}$ (масштабування).
  3. $E_{i \leftarrow R_i + c R_j}^{-1} = E_{i \leftarrow R_i - c R_j}$ (елімінація).

### 2.2. Оборотні та необоротні матриці (Invertible / Non-singular Matrices)
Квадратна матриця $A \in \mathbb{R}^{n \times n}$ називається **оборотною (invertible / non-singular)**, якщо існує матриця $A^{-1} \in \mathbb{R}^{n \times n}$ така, що:
$$A A^{-1} = A^{-1} A = I_n$$

#### Теорема про оборотну матрицю (Invertible Matrix Theorem):
Для квадратної матриці $A \in \mathbb{R}^{n \times n}$ наступні твердження є **еквівалентними**:
1. $A$ є оборотною ($A^{-1}$ існує).
2. $\text{rank}(A) = n$ (повний ранг / full rank).
3. $\text{RREF}(A) = I_n$.
4. $\det(A) \neq 0$.
5. Рівняння $A\mathbf{x} = \mathbf{0}$ має лише тривіальний розв'язок $\mathbf{x} = \mathbf{0}$ ($\text{null}(A) = \{\mathbf{0}\}$).
6. Рівняння $A\mathbf{x} = \mathbf{b}$ має єдиний розв'язок для будь-якого $\mathbf{b} \in \mathbb{R}^n$: $\mathbf{x} = A^{-1}\mathbf{b}$.
7. Стовпці (і рядки) матриці $A$ є лінійно незалежними.
8. Жодне з власних значень $\lambda_i \neq 0$.

### 2.3. Алгоритм Гаусса-Жордана знаходження $A^{-1}$
Формуємо блочну матрицю $[A \mid I_n]$ і зводимо її до RREF:
$$[A \mid I_n] \xrightarrow{\text{Gauss-Jordan}} [I_n \mid A^{-1}]$$
Якщо зліва не вдається отримати $I_n$ (з'явився нульовий рядок), матриця $A$ є **сингулярною (необоротною)**.

---

## 3. LU та LUP факторизація (LU & LUP Factorization)

### 3.1. Суть LU-розкладу
Будь-яку квадратну матрицю $A$ (якщо процес виключення Гаусса не вимагає перестановки рядків) можна розкласти у добуток:
$$A = L U$$
де:
* $L \in \mathbb{R}^{n \times n}$ — **нижня трикутна матриця (Lower triangular matrix)** з одиницями на головній діагоналі (Unit Lower Triangular):
  $$L = \begin{bmatrix} 1 & 0 & \dots & 0 \\ l_{21} & 1 & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ l_{n1} & l_{n2} & \dots & 1 \end{bmatrix}$$
  де $l_{ik} = \frac{a_{ik}^{(k-1)}}{a_{kk}^{(k-1)}}$ — множники (multipliers), що використовувались для обнулення елементів під час прямого ходу Гаусса.
* $U \in \mathbb{R}^{n \times n}$ — **верхня трикутна матриця (Upper triangular matrix)**, отримана в результаті прямого ходу виключення Гаусса (REF матриці $A$):
  $$U = \begin{bmatrix} u_{11} & u_{12} & \dots & u_{1n} \\ 0 & u_{22} & \dots & u_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & u_{nn} \end{bmatrix}$$

### 3.2. Розв'язання $A\mathbf{x} = \mathbf{b}$ за допомогою LU
Замість $A\mathbf{x} = \mathbf{b}$ розв'язуємо дві прості трикутні системи:
1. **Пряма підстановка (Forward substitution):** $L \mathbf{y} = \mathbf{b}$ (складність $\mathcal{O}(n^2)$).
2. **Зворотна підстановка (Back substitution):** $U \mathbf{x} = \mathbf{y}$ (складність $\mathcal{O}(n^2)$).

> **💡 Чому це вигідно в інженерії та ML?**
> Знаходження LU-розкладу займає $\mathcal{O}(\frac{2}{3} n^3)$ операцій. Якщо нам потрібно розв'язати сотні систем $A\mathbf{x} = \mathbf{b}_k$ з тією самою матрицею $A$ і різними правими частинами $\mathbf{b}_k$, ми робимо розклад $A=LU$ **один раз**, а потім кожен новий вектор розв'язуємо за швидкі $\mathcal{O}(n^2)$ операцій замість повторного виключення Гаусса!

### 3.3. LUP-розклад (з перестановкою рядків / Partial Pivoting)
Якщо під час виключення на діагоналі з'являється 0 (або для числової стійкості від переповнення), рядки переставляють за допомогою **матриці перестановки (Permutation matrix $P$)**:
$$P A = L U$$
де $P$ — матриця, отримана перестановкою рядків з $I_n$ ($P^{-1} = P^T$).

---

## 4. ⚠️ Типові підводні камені (Common Pitfalls)
* 🚨 **$(AB)^{-1} = B^{-1} A^{-1}$ (порядок змінюється!):** Аналогічно до $(AB)^T = B^T A^T$.
* 🚨 **Неіснування $A^{-1}$ для прямокутних матриць:** Звичайна двостороння обернена матриця $A^{-1}$ визначена **лише для квадратних** матриць повного рангу. Для прямокутних матриць використовують псевдообернену матрицю Мура-Пенроуза ($A^+$).

---

## 📊 Зведена таблиця (Summary Table)

| Операція / Поняття | Формула | Властивості |
|---|---|---|
| Множення матриць | $C = AB, \quad c_{ij} = \sum a_{ik} b_{kj}$ | Асоціативне, але $AB \neq BA$ |
| Транспонування | $(AB)^T = B^T A^T$ | $(A^T)^{-1} = (A^{-1})^T$ |
| Обернена матриця | $A A^{-1} = I$ | $(AB)^{-1} = B^{-1} A^{-1}$ |
| LU-факторизація | $A = LU$ | $L$ — lower unit triangular, $U$ — upper triangular |
| LUP-факторизація | $PA = LU$ | $P$ — permutation matrix ($P^T P = I$) |

---

[⬅️ 01. Системи лінійних рівнянь](01-systems-of-linear-equations-and-gaussian-elimination.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 03. Визначники ➡️](03-determinants-and-cramers-rule.md) | [⚡ Cheat Sheet](cheat-sheet.md)
