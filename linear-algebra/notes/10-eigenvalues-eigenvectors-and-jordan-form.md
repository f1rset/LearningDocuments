# 10. Власні значення, власні вектори, діагоналізація та жорданова форма (Eigenvalues, Eigenvectors, Diagonalization & Jordan Normal Form)

[⬅️ 09. Ортогональні матриці та QR-розклад](09-orthogonal-matrices-gram-schmidt-and-qr.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 11. Симетричні матриці та квадратичні форми ➡️](11-symmetric-matrices-quadratic-forms-and-pca.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Власні значення та власні вектори (Eigenvalues & Eigenvectors):** Напрямки, в яких лінійне перетворення діє як просте розтягування/стискання: $A\mathbf{v} = \lambda \mathbf{v}$.
2. **Характеристичний многочлен (Characteristic Polynomial):** $\det(A - \lambda I) = 0$.
3. **Діагоналізація матриць (Matrix Diagonalization):** Розклад $A = S \Lambda S^{-1}$.
4. **Алгебраїчна та геометрична кратність (Algebraic & Geometric Multiplicity):** Критерій діагоналізовності.
5. **Жорданова нормальна форма (Jordan Normal Form / JNF):** Канонічний вигляд недіагоналізовних матриць.
6. **Застосування:** Різницеві рівняння (динамічні системи, ряди Маркова) та системи лінійних диференціальних рівнянь.

---

## 1. Власні значення та вектори (Eigenvalues & Eigenvectors)

### 1.1. Означення
Нехай $A \in \mathbb{R}^{n \times n}$ — квадратна матриця.
Ненульовий вектор $\mathbf{v} \neq \mathbf{0}$ називається **власним вектором (eigenvector)** матриці $A$, що відповідає **власному значенню $\lambda$ (eigenvalue)** ($\lambda \in \mathbb{C}$), якщо:
$$\mathbf{A \mathbf{v} = \lambda \mathbf{v} \iff (A - \lambda I_n) \mathbf{v} = \mathbf{0}}$$

> **💡 Геометрична інтуїція (Geometric Intuition):**
> Зазвичай дія матриці $A$ на довільний вектор $\mathbf{x}$ повертає його і змінює його довжину. Але для **власного вектора $\mathbf{v}$** матриця діє як **звичайний скаляр**: вона **НЕ змінює його напрямок**, а лише масштабує його у $\lambda$ разів!
> * Якщо $\lambda > 1$ — вектор розтягується.
> * Якщо $0 < \lambda < 1$ — стискається.
> * Якщо $\lambda < 0$ — напрямок змінюється на протилежний.
> * Якщо $\lambda = 0$ — вектор стискається в нуль ($\mathbf{v} \in \text{Null}(A)$).

```text
       Довільний вектор x                   Власний вектор v (A*v = lambda*v)
              ^                                    ^
          A*x |   / x                          A*v |
              |  / (Поворот + розтяг)              |  / v (Тільки розтяг вздовж
              | /                                  | /     тієї самої прямої!)
             -+------>                            -+------>
```

### 1.2. Алгоритм знаходження власних значень і векторів

1. **Характеристичне рівняння (Characteristic Equation):**
   Оскільки $\mathbf{v} \neq \mathbf{0}$, система $(A - \lambda I)\mathbf{v} = \mathbf{0}$ має нетривіальні розв'язки $\iff$ матриця $(A - \lambda I)$ є виродженою:
   $$p_A(\lambda) = \det(A - \lambda I_n) = 0$$
   * $p_A(\lambda)$ — многочлен степеня $n$. Його $n$ коренів $\lambda_1, \lambda_2, \dots, \lambda_n \in \mathbb{C}$ є власними значеннями $A$.
2. **Власний підпростір (Eigenspace $E_\lambda$):**
   Для кожного знайденого $\lambda_i$ власний вектор знаходиться як нетривіальний розв'язок однорідної системи:
   $$E_{\lambda_i} = \text{Null}(A - \lambda_i I) = \{\mathbf{v} \in \mathbb{R}^n \mid (A - \lambda_i I)\mathbf{v} = \mathbf{0}\}$$

### 1.3. Властивості власних значень:
* **Слід матриці (Trace):** Сума всіх власних значень дорівнює сліду:
  $$\sum_{i=1}^n \lambda_i = \text{tr}(A) = \sum_{i=1}^n a_{ii}$$
* **Визначник (Determinant):** Добуток усіх власних значень дорівнює визначнику:
  $$\prod_{i=1}^n \lambda_i = \det(A)$$
* Матриця $A$ **оборотна** $\iff$ усі $\lambda_i \neq 0$.
* Власні значення $A^k$ дорівнюють $\lambda_i^k$, власні вектори залишаються тими самими $\mathbf{v}_i$.
* Власні значення $A^{-1}$ дорівнюють $1/\lambda_i$.

---

## 2. Діагоналізація матриць (Matrix Diagonalization)

### 2.1. Теорема про діагоналізацію
Квадратна матриця $A \in \mathbb{R}^{n \times n}$ називається **діагоналізовною (diagonalizable)**, якщо вона подібна до діагональної матриці $\Lambda$:
$$\mathbf{A = S \Lambda S^{-1} \iff S^{-1} A S = \Lambda}$$
де:
* $S = [\mathbf{v}_1 \mid \mathbf{v}_2 \mid \dots \mid \mathbf{v}_n] \in \mathbb{R}^{n \times n}$ — **матриця власних векторів (eigenvector matrix)**.
* $\Lambda = \text{diag}(\lambda_1, \lambda_2, \dots, \lambda_n)$ — **діагональна матриця власних значень (eigenvalue matrix)**.

$$
A \begin{bmatrix} \mathbf{v}_1 & \dots & \mathbf{v}_n \end{bmatrix} = \begin{bmatrix} \lambda_1 \mathbf{v}_1 & \dots & \lambda_n \mathbf{v}_n \end{bmatrix} = \begin{bmatrix} \mathbf{v}_1 & \dots & \mathbf{v}_n \end{bmatrix} \begin{bmatrix} \lambda_1 & & 0 \\ & \ddots & \\ 0 & & \lambda_n \end{bmatrix} \implies AS = S\Lambda
$$

### 2.2. Критерій діагоналізовності
Матриця $A \in \mathbb{R}^{n \times n}$ є діагоналізовною тоді і тільки тоді, коли вона має **$n$ лінійно незалежних власних векторів**.

* **Достатня умова:** Якщо всі $n$ власних значень $\lambda_1, \dots, \lambda_n$ є **різними (distinct)**, матриця $A$ **гарантовано діагоналізовна**.
* **Випадок кратних коренів:**
  * **Алгебраїчна кратність (Algebraic Multiplicity / $AM(\lambda)$):** Кратність кореня $\lambda$ у характеристичному многочлені $p_A(\lambda)$.
  * **Геометрична кратність (Geometric Multiplicity / $GM(\lambda)$):** Кількість лінійно незалежних власних векторів для даного $\lambda$: $GM(\lambda) = \dim(\text{Null}(A - \lambda I))$.
  * Завжди: $1 \le GM(\lambda) \le AM(\lambda)$.
  * **Критерій:** Матриця діагоналізовна $\iff GM(\lambda_i) = AM(\lambda_i)$ для всіх власних значень.

---

## 3. Обчислення степенів матриці (Powers of a Matrix)

Якщо матриця $A$ діагоналізовна ($A = S\Lambda S^{-1}$), то піднесення до степеня $k$ стає тривіальним:
$$A^k = (S \Lambda S^{-1})(S \Lambda S^{-1}) \dots (S \Lambda S^{-1}) = \mathbf{S \Lambda^k S^{-1}}$$
$$
\Lambda^k = \begin{bmatrix} \lambda_1^k & & 0 \\ & \ddots & \\ 0 & & \lambda_n^k \end{bmatrix}
$$
*(Складність обчислення $A^{1000}$ знижується з $1000 \cdot \mathcal{O}(n^3)$ до одного розкладу $\mathcal{O}(n^3)$!).*

---

## 4. Жорданова нормальна форма (Jordan Normal Form / JNF)

### 4.1. Що робити, якщо матриця недіагоналізовна ($GM < AM$)?
Якщо власних векторів не вистачає для утворення базису (дефектна матриця / defective matrix), матрицю не можна звести до чисто діагонального вигляду. Проте її завжди можна звести над полем $\mathbb{C}$ до майже діагональної **Жорданової форми**:
$$A = M J M^{-1}$$
де $J$ — блочно-діагональна матриця:
$$J = \begin{bmatrix} J_1 & & 0 \\ & \ddots & \\ 0 & & J_k \end{bmatrix}, \quad J_i(\lambda) = \begin{bmatrix} \lambda & 1 & 0 & \dots & 0 \\ 0 & \lambda & 1 & \dots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \dots & \lambda & 1 \\ 0 & 0 & \dots & 0 & \lambda \end{bmatrix}$$
* $J_i(\lambda)$ — **жорданова клітина (Jordan block)** розміру $m \times m$ з власним значенням $\lambda$ на діагоналі та **одиницями над головною діагоналлю (superdiagonal)**.
* Кількість жорданових клітин для значення $\lambda$ в точності дорівнює його геометричній кратності $GM(\lambda)$.
* Для побудови матриці переходу $M$ використовують **приєднані (узагальнені) власні вектори (generalized eigenvectors)**: $(A - \lambda I)\mathbf{w}_2 = \mathbf{v}_1$.

---

## 5. Застосування до динамічних систем та диференціальних рівнянь

### 5.1. Лінійні різницеві рівняння (Discrete Dynamical Systems / Markov Chains)
$\mathbf{x}_{k+1} = A \mathbf{x}_k \implies \mathbf{x}_k = A^k \mathbf{x}_0 = S \Lambda^k S^{-1} \mathbf{x}_0 = \sum_{i=1}^n c_i \lambda_i^k \mathbf{v}_i$.
* **Асимптотична стійкість (Stability):**
  * Якщо всі $|\lambda_i| < 1 \implies \mathbf{x}_k \to \mathbf{0}$ при $k \to \infty$ (стійкий фокус / атрактор).
  * Якщо хоча б один $|\lambda_i| > 1 \implies$ система експоненційно розбігається ($\mathbf{x}_k \to \infty$).
  * Якщо $\lambda_1 = 1$, а решта $|\lambda_i| < 1 \implies$ система збігається до стаціонарного стану (Steady State / Марковські ланцюги, Google PageRank).

### 5.2. Системи лінійних диференціальних рівнянь (System of ODEs)
$$\frac{d\mathbf{x}}{dt} = A \mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0$$
* **Матрична експонента (Matrix Exponential $e^{At}$):**
  $$\mathbf{x}(t) = e^{At} \mathbf{x}_0 = S e^{\Lambda t} S^{-1} \mathbf{x}_0 = \sum_{i=1}^n c_i e^{\lambda_i t} \mathbf{v}_i$$
  де $e^{\Lambda t} = \text{diag}(e^{\lambda_1 t}, e^{\lambda_2 t}, \dots, e^{\lambda_n t})$.
* **Стійкість неперервної системи:** Система стійка $\iff \text{Re}(\lambda_i) < 0$ для всіх власних значень.

---

## 📊 Зведена таблиця (Summary Table)

| Поняття | Формула | Умова / Властивість |
|---|---|---|
| Власне рівняння | $A\mathbf{v} = \lambda \mathbf{v}$ | $\mathbf{v} \neq \mathbf{0}$ |
| Характеристичний многочлен | $\det(A - \lambda I) = 0$ | Степінь $n$, дає $n$ власних значень |
| Діагоналізація | $A = S\Lambda S^{-1}$ | Потрібно $n$ лінійно незалежних власних векторів |
| Степінь матриці | $A^k = S\Lambda^k S^{-1}$ | Швидке піднесення до степеня |
| Жорданова форма | $A = MJM^{-1}$ | Існує завжди для будь-якої квадратної матриці |
| Диференціальне рівняння | $\mathbf{x}(t) = \sum c_i e^{\lambda_i t} \mathbf{v}_i$ | Розв'язок $\dot{\mathbf{x}} = A\mathbf{x}$ |

---

[⬅️ 09. Ортогональні матриці та QR-розклад](09-orthogonal-matrices-gram-schmidt-and-qr.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 11. Симетричні матриці та квадратичні форми ➡️](11-symmetric-matrices-quadratic-forms-and-pca.md) | [⚡ Cheat Sheet](cheat-sheet.md)

