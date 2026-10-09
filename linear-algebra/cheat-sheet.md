# ⚡ Ultimate Linear Algebra Cheat Sheet: Вся лінійна алгебра в формулах та концепціях

[🏠 Головний зміст](index.md) | [01. Системи лінійних рівнянь](01-systems-of-linear-equations-and-gaussian-elimination.md) | [13. Сингулярний розклад (SVD)](13-svd-and-best-rank-approximations.md)

---

## 🔑 1. Фундаментальні концепції та матричні розклади (Matrix Factorizations Overview)

```text
               Шпаргалка головних факторизацій матриць (Matrix Factorizations):
               
   Розклад       Формула             Вимоги до матриці        Сфери застосування
   -----------------------------------------------------------------------------------------
   LU            A = L U             Квадратна                Швидке розв'язання Ax = b
   LUP           P A = L U           Квадратна                Чисельно стійке виключення Гаусса
   Cholesky      A = L L^T           Симетрична A = A^T > 0   МНК, оптимізація, фільтри Калмана
   QR            A = Q R             Прямокутна (m >= n)      Чисельно стійкий МНК, ортогоналізація
   Spectral      A = Q Λ Q^T         Симетрична A = A^T       PCA, аналіз квадратичних форм
   Eigendecomp   A = S Λ S^-1        n лінійно незал. векторів Динамічні системи, матрична експонента e^(At)
   Jordan        A = M J M^-1        Будь-яка квадратна       Аналіз недіагоналізовних систем
   SVD           A = U Σ V^T         БУДЬ-ЯКА матриця m x n   Стиснення, псевдообернення, рекомендації
```

---

## ⚡ 2. Зведена таблиця формул та алгоритмів

| Розділ / Тема | Ключова формула / Співвідношення | Англійський термін | Пояснення / Фізичний зміст |
|---|---|---|---|
| **Системи рівнянь** | $A\mathbf{x} = \mathbf{b}$ | System of Linear Equations | $m$ рівнянь, $n$ невідомих |
| **Теорема Кронекера-Капеллі** | $\text{rank}(A) = \text{rank}([A \mid \mathbf{b}])$ | Rouché–Capelli Theorem | Критерій сумісності системи |
| **Добуток матриць** | $c_{ij} = \sum_{k} a_{ik} b_{kj}$ | Matrix Multiplication | Асоціативне, але $AB \neq BA$ |
| **Транспонування** | $(AB)^T = B^T A^T$ | Transpose | Порядок співмножників змінюється |
| **Обернена матриця** | $(AB)^{-1} = B^{-1} A^{-1}$ | Inverse Matrix | Визначена для $\det(A) \neq 0$ |
| **Визначник $2 \times 2$** | $\det \begin{bmatrix} a & b \\ c & d \end{bmatrix} = ad - bc$ | Determinant | Орієнтована площа |
| **Властивість визначника** | $\det(AB) = \det(A)\det(B)$ | Multiplicative Property | $\det(A^{-1}) = 1/\det(A)$ |
| **Правило Крамера** | $x_i = \frac{\det(A_i)}{\det(A)}$ | Cramer's Rule | Розв'язок через замінені стовпці |
| **Ранг і дефект** | $\text{rank}(A) + \dim(\text{Null}(A)) = n$ | Rank-Nullity Theorem | Фундаментальна теорема |
| **4 підпростори** | $\text{Null}(A) = (\text{Row}(A))^\perp$ | Orthogonal Complements | $\text{Null}(A^T) = (\text{Col}(A))^\perp$ |
| **Зміна базису** | $B = P^{-1} A P$ | Similar Matrices | $P$ — матриця переходу |
| **Скалярний добуток** | $\mathbf{u}^T \mathbf{v} = \|\mathbf{u}\|_2 \|\mathbf{v}\|_2 \cos\theta$ | Dot / Inner Product | Міра кута та проекції |
| **Коші-Буняковський** | $|\mathbf{u}^T \mathbf{v}| \le \|\mathbf{u}\|_2 \|\mathbf{v}\|_2$ | Cauchy-Schwarz Inequality | Рівність $\iff \mathbf{u} = c\mathbf{v}$ |
| **Проекція на підпростір** | $P = A (A^T A)^{-1} A^T$ | Projection Matrix | $P^2 = P, \quad P^T = P$ |
| **Нормальні рівняння** | $A^T A \hat{\mathbf{x}} = A^T \mathbf{b}$ | Normal Equations | Метод найменших квадратів (OLS) |
| **МНК через QR** | $R \hat{\mathbf{x}} = Q^T \mathbf{b}$ | Least Squares via QR | Чисельно стійкий алгоритм |
| **Власні значення** | $A\mathbf{v} = \lambda \mathbf{v}, \quad \det(A - \lambda I) = 0$ | Eigenvalues & Eigenvectors | $\sum \lambda_i = \text{tr}(A), \prod \lambda_i = \det(A)$ |
| **Діагоналізація** | $A = S \Lambda S^{-1}$ | Diagonalization | $S$ — матриця власних векторів |
| **Спектральна теорема** | $A = Q \Lambda Q^T = \sum \lambda_i \mathbf{q}_i \mathbf{q}_i^T$ | Spectral Theorem (Symmetric) | Усі $\lambda_i \in \mathbb{R}$, $Q^T Q = I$ |
| **Критерій Сильвестра** | $\Delta_1 > 0, \Delta_2 > 0, \dots, \Delta_n > 0$ | Sylvester's Criterion | Додатна визначеність $A \succ 0$ |
| **Розклад Холецького** | $A = L L^T$ | Cholesky Decomposition | $L$ — lower triangular ($l_{ii} > 0$) |
| **Сингулярний розклад (SVD)** | $A = U \Sigma V^T = \sum \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ | SVD | $\sigma_i = \sqrt{\lambda_i(A^T A)}$ |
| **Псевдообернена Мура-Пенроуза** | $A^+ = V \Sigma^+ U^T$ | Pseudoinverse | $\hat{\mathbf{x}} = A^+ \mathbf{b}$ |
| **Теорема Екхарта-Янга** | $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ | Eckart-Young-Mirsky Theorem | Найкраще наближення рангу $k$ |

---

## 🛠️ 3. Інженерні та математичні правила (Rules of Thumb)

* 🔹 **Повний стовпчиковий ранг:** Якщо матриця $A \in \mathbb{R}^{m \times n}$ має $\text{rank}(A) = n$ ($m \ge n$), то $A^T A$ є квадратною строго **додатно визначеною ($A^TA \succ 0$)** і завжди **оборотною**.
* 🔹 **Діагоналізація vs SVD:** Власні значення $\lambda_i$ показують поведінку оператора при багаторазовому застосуванні ($A^k$, динамічні системи). Сингулярні числа $\sigma_i$ показують геометричне розтягнення осей та енергію матриці для прямокутних даних (стиснення, PCA).
* 🔹 **Чисельна стійкість:** Для розв'язання $A\mathbf{x} = \mathbf{b}$ **ніколи не рахуйте $A^{-1}$ явно через визначник**! Використовуйте $LU$-розклад для квадратних систем і $QR$-розклад для перевизначених систем найменших квадратів.

---

## ⚠️ 4. Критичні підводні камені (Common Pitfalls)

* 🚨 **Некомутативність матриць:** $AB \neq BA$. Відповідно, $(A+B)^2 \neq A^2 + 2AB + B^2$ (правильно: $A^2 + AB + BA + B^2$).
* 🚨 **Власні значення суми чи добутку:** У загальному випадку $\lambda(A + B) \neq \lambda(A) + \lambda(B)$ та $\lambda(AB) \neq \lambda(A)\lambda(B)$ (це виконується лише якщо $A$ та $B$ комутують, тобто мають спільний базис власних векторів).
* 🚨 **Простори стовпців і рядків:** Рядкові перетворення Гаусса зберігають простір рядків $\text{Row}(A)$ та ядро $\text{Null}(A)$, але **змінюють простір стовпців $\text{Col}(A)$**. Базис $\text{Col}(A)$ завжди береться з початкової матриці $A$!

---

[🏠 Назад до Змісту курсу](index.md)

