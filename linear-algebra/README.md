# 📐 Лінійна алгебра: Повний курс-конспект (Linear Algebra Complete Course)

[🏠 Головна сторінка репозиторію](../../README.md) | [⚡ Ultimate Cheat Sheet](cheat-sheet.md)

Цей курс містить вичерпні, глибокі та структуровані матеріали з **Лінійної алгебри та її застосувань у Data Science, Machine Learning та інженерії**. Всі терміни та назви теорем повністю дублюються **англійською мовою**.

---

> **🔥 Швидкий доступ:** Повна інженерна шпаргалка з усіма матричними розкладами ($LU, QR, SVD, \text{Cholesky}$), формулами та теоремами доступна тут 👉 [**⚡ ULTIMATE LINEAR ALGEBRA CHEAT SHEET**](cheat-sheet.md)

---

## 🗺️ Навігація по модулях курсу (Modules & Topics)

| № | Тема / Модуль | Ключові концепції (Core Concepts) | Англійські терміни | Посилання |
|---|---|---|---|---|
| ⭐ | **Cheat Sheet** | **Усі факторизації ($LU, QR, SVD$), теореми, формули та підводні камені** | **Ultimate Cheat Sheet** | [**⚡ Відкрити Cheat Sheet**](cheat-sheet.md) |
| **01** | **Системи лінійних рівнянь та метод Гаусса** | Метод виключення Гаусса та Гаусса-Жордана, REF, RREF, теорема Кронекера-Капеллі, існування та єдиність. | *Systems of Linear Equations, Gaussian & Gauss-Jordan Elimination, REF, RREF, Rouché–Capelli Theorem* | [Читати модуль 01 ➡️](notes/01-systems-of-linear-equations-and-gaussian-elimination.md) |
| **02** | **Матриці, вектори та LU-розклад** | Базові операції, елементарні матриці, критерії оборотності (Invertible Matrix Theorem), $LU$ та $LUP$ факторизація. | *Matrices & Vectors, Elementary Transformations, Invertible Matrices, LU & LUP Factorization* | [Читати модуль 02 ➡️](notes/02-matrices-vectors-and-lu-factorization.md) |
| **03** | **Визначники та правило Крамера** | Визначники 2D/3D, розклад Лапласа, властивості $\det$, орієнтований об'єм, векторний добуток, союзна матриця, правило Крамера. | *Determinants in 2D & 3D, Laplace Expansion, Properties, Volume & Cross Product, Cramer's Rule* | [Читати модуль 03 ➡️](notes/03-determinants-and-cramers-rule.md) |
| **04** | **Векторні простори та підпростори** | Аксіоми лінійного простору, критерій підпростору, лінійні комбінації, лінійна оболонка ($\text{span}$), лінійна незалежність. | *Vector Spaces, Subspaces, Linear Combinations, Span, Linear Independence* | [Читати модуль 04 ➡️](notes/04-vector-spaces-subspaces-and-spans.md) |
| **05** | **Базиси та 4 фундаментальні підпростори** | Базис, розмірність, ізоморфізм $\mathbb{R}^n$, $\text{Col}(A), \text{Row}(A), \text{Null}(A), \text{Null}(A^T)$, теорема про ранг і дефект. | *Bases, Dimension, Isomorphism, The Four Fundamental Subspaces, Rank-Nullity Theorem* | [Читати модуль 05 ➡️](notes/05-bases-dimension-rank-and-four-subspaces.md) |
| **06** | **Заміна базису та лінійні оператори** | Матриця переходу $P$, матриця лінійного відображення, подібні матриці ($B = P^{-1} A P$), інваріанти оператора. | *Change of Basis, Transition Matrices, Linear Transformations, Similar Matrices* | [Читати модуль 06 ➡️](notes/06-change-of-basis-and-linear-transformations.md) |
| **07** | **Норми, скалярний добуток та ортогональність** | Аксіоми норми ($L_1, L_2, L_\infty, L_p$), нерівність Коші-Буняковського-Шварца, кут між векторами, ортогональне доповнення $W^\perp$. | *Norms, Inner Products, Cauchy-Schwarz Inequality, Angles, Orthogonal Complements* | [Читати модуль 07 ➡️](notes/07-norms-inner-products-and-orthogonality.md) |
| **08** | **Ортонормовані базиси, проекції та МНК** | Ортонормовані базиси (ONB), матриця проекції $P = A(A^TA)^{-1}A^T$, нормальні рівняння, метод найменших квадратів (OLS). | *Orthonormal Bases, Orthogonal Projections, Projection Matrix, Normal Equations, Least Squares* | [Читати модуль 08 ➡️](notes/08-orthonormal-bases-projections-and-least-squares.md) |
| **09** | **Ортогональні матриці, Грама-Шмідта та QR-розклад** | Ортогональні матриці ($Q^T Q = I$), ізометрія, алгоритм Грама-Шмідта (CGS/MGS), $A = QR$, чисельно стійкий МНК ($R\hat{\mathbf{x}} = Q^T\mathbf{b}$). | *Orthogonal Matrices, Gram-Schmidt Process, QR Factorization, Stable Least Squares* | [Читати модуль 09 ➡️](notes/09-orthogonal-matrices-gram-schmidt-and-qr.md) |
| **10** | **Власні значення, вектори та жорданова форма** | $A\mathbf{v} = \lambda \mathbf{v}$, характеристичне рівняння, діагоналізація $A = S\Lambda S^{-1}$, жорданова нормальна форма (JNF), динамічні системи. | *Eigenvalues & Eigenvectors, Diagonalization, Jordan Normal Form, Dynamical Systems, Matrix Exponential* | [Читати модуль 10 ➡️](notes/10-eigenvalues-eigenvectors-and-jordan-form.md) |
| **11** | **Симетричні матриці, квадратичні форми, Холецький та PCA** | Спектральна теорема ($A = Q\Lambda Q^T$), додатна визначеність ($A \succ 0$), критерій Сильвестра, розклад Холецького ($A = LL^T$), метод PCA. | *Symmetric Matrices, Quadratic Forms, Positive Definiteness, Cholesky Decomposition, PCA* | [Читати модуль 11 ➡️](notes/11-symmetric-matrices-quadratic-forms-and-pca.md) |
| **12** | **Унітарні матриці та спектральна теорія** | Кососиметричні, ермітові, унітарні та нормальні матриці ($A^* A = A A^*$), теорема Шура, розклад одиниці через проектори ($I = \sum P_i$). | *Skew-Symmetric, Hermitian, Unitary & Normal Matrices, Schur Theorem, Spectral Decomposition* | [Читати модуль 12 ➡️](notes/12-symmetric-unitary-matrices-and-spectral-decomposition.md) |
| **13** | **Сингулярний розклад (SVD) та застосування** | $A = U\Sigma V^T$, сингулярні числа $\sigma_i$, теорема Екхарта-Янга (Low-Rank), стиснення зображень, псевдообернена Мура-Пенроуза ($A^+$). | *Singular Value Decomposition (SVD), Low-Rank Approximation, Eckart-Young Theorem, Pseudoinverse, Applications* | [Читати модуль 13 ➡️](notes/13-svd-and-best-rank-approximations.md) |

---

[🏠 Головна сторінка репозиторію](../../README.md)
