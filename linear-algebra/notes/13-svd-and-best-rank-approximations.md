# 13. Сингулярний розклад, наближення низького рангу та застосування SVD (Singular Value Decomposition & Applications)

[⬅️ 12. Спектральна теорія та унітарні матриці](12-symmetric-unitary-matrices-and-spectral-decomposition.md) | [🏠 Головний зміст](index.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Сингулярний розклад (Singular Value Decomposition / SVD):** Вершина лінійної алгебри: $A = U \Sigma V^T$ для **будь-якої** прямокутної матриці $m \times n$.
2. **Сингулярні числа (Singular Values $\sigma_i$):** $\sigma_i = \sqrt{\lambda_i(A^T A)} \ge 0$.
3. **Лівий та правий сингулярні вектори (Left & Right Singular Vectors $U$ and $V$):** Ортонормовані базиси 4 фундаментальних підпросторів.
4. **Теорема Екхарта-Янга-Мірського (Eckart-Young-Mirsky Theorem):** Найкраще наближення матриці низького рангу $k$ (Low-Rank Approximation).
5. **Прикладні застосування SVD:** Стиснення зображень, рекомендаційні системи (Matrix Factorization), псевдообернення Мура-Пенроуза ($A^+$), латентно-семантичний аналіз (LSA / LSI).

---

## 1. Сингулярний розклад (Singular Value Decomposition / SVD)

### 1.1. Головна теорема SVD
Для **будь-якої** дійсної прямокутної матриці $A \in \mathbb{R}^{m \times n}$ рангу $\text{rank}(A) = r \le \min(m, n)$ існує факторизація:
$$\mathbf{A = U \Sigma V^T = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T}$$
де:
* $U \in \mathbb{R}^{m \times m}$ — **ортогональна матриця лівих сингулярних векторів (Left Singular Vectors)** ($U^T U = I_m$).
* $V \in \mathbb{R}^{n \times n}$ — **ортогональна матриця правих сингулярних векторів (Right Singular Vectors)** ($V^T V = I_n$).
* $\Sigma \in \mathbb{R}^{m \times n}$ — **прямокутна діагональна матриця сингулярних чисел (Singular Values)**:
  $$\Sigma = \begin{bmatrix} \sigma_1 & & & 0 \\ & \sigma_2 & & 0 \\ & & \ddots & \vdots \\ 0 & 0 & \dots & \sigma_r \\ 0 & 0 & \dots & 0 \end{bmatrix}, \quad \sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_r > 0$$

```text
       Геометрична інтуїція SVD (Дія матриці як 3 етапи):
       
       x в R^n              V^T x                 Sigma V^T x            A x = U Sigma V^T x в R^m
        (Одинична сфера)     (Обертання)           (Розтягнення осей)     (Обертання в R^m)
            .-.                   .-.                    .---.                   .---.
           (   )     ===>        (   )     ===>         (  *  )   ===>          (  /  )
            `-'                   `-'                    `---'                   `---'
                          Обертання базису      Півосі еліпсоїда       Остаточний еліпсоїд
                          правими векторами v_i  масштабуються на σ_i  орієнтований векторами u_i
```

---

## 2. Як знаходяться $U$, $\Sigma$ та $V$?

### 2.1. Зв'язок із симетричними матрицями $A^T A$ та $A A^T$
Розглянемо добутки матриці $A$ на її транспоновану:
1. **$A^T A \in \mathbb{R}^{n \times n}$ (Симетрична, додатно напіввизначена):**
   $$A^T A = (U \Sigma V^T)^T (U \Sigma V^T) = V \Sigma^T \underbrace{U^T U}_{I_m} \Sigma V^T = \mathbf{V (\Sigma^T \Sigma) V^T}$$
   * Стовпці матриці $V$ — це **власні вектори матриці $A^T A$**.
   * Сингулярні числа — це квадратні корені з ненульових власних значень $A^T A$:
     $$\mathbf{\sigma_i = \sqrt{\lambda_i(A^T A)}}$$
2. **$A A^T \in \mathbb{R}^{m \times m}$ (Симетрична, додатно напіввизначена):**
   $$A A^T = (U \Sigma V^T)(U \Sigma V^T)^T = U \Sigma \underbrace{V^T V}_{I_n} \Sigma^T U^T = \mathbf{U (\Sigma \Sigma^T) U^T}$$
   * Стовпці матриці $U$ — це **власні вектори матриці $A A^T$**.
3. **Зв'язок між векторами:**
   $$A \mathbf{v}_i = \sigma_i \mathbf{u}_i \implies \mathbf{u}_i = \frac{1}{\sigma_i} A \mathbf{v}_i$$

### 2.2. SVD та 4 фундаментальні підпростори
SVD дає **ідеальні ортонормовані базиси для всіх 4 підпросторів матриці $A$**:

| Підпростір (Subspace) | Розмірність | Ортонормований базис |
|---|---|---|
| **Простір рядків $\text{Row}(A)$** | $r$ | Перші $r$ стовпців матриці $V$: $\{\mathbf{v}_1, \dots, \mathbf{v}_r\}$ |
| **Ядро $\text{Null}(A)$** | $n - r$ | Останні $n-r$ стовпців матриці $V$: $\{\mathbf{v}_{r+1}, \dots, \mathbf{v}_n\}$ |
| **Простір стовпців $\text{Col}(A)$** | $r$ | Перші $r$ стовпців матриці $U$: $\{\mathbf{u}_1, \dots, \mathbf{u}_r\}$ |
| **Ліве ядро $\text{Null}(A^T)$** | $m - r$ | Останні $m-r$ стовпців матриці $U$: $\{\mathbf{u}_{r+1}, \dots, \mathbf{u}_m\}$ |

---

## 3. Наближення матриці низького рангу (Low-Rank Approximation)

### 3.1. Розклад SVD як сума матриць рангу 1 (Outer Product Form)
$$A = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T = \sigma_1 \mathbf{u}_1 \mathbf{v}_1^T + \sigma_2 \mathbf{u}_2 \mathbf{v}_2^T + \dots + \sigma_r \mathbf{u}_r \mathbf{v}_r^T$$
Кожен доданок $\mathbf{u}_i \mathbf{v}_i^T \in \mathbb{R}^{m \times n}$ є матрицею **рангу 1**.

### 3.2. Теорема Екхарта-Янга-Мірського (Eckart-Young-Mirsky Theorem)
Нехай ми хочемо наблизити матрицю $A$ рангу $r$ іншою матрицею $A_k$ значно меншого рангу $k < r$ ($k$-рангова апроксимація).

**Найкращою матрицею рангу $k$** у сенсі як спектральної норми $\| \cdot \|_2$, так і норми Фробеніуса $\| \cdot \|_F$ є **усічена сума SVD (Truncated SVD)**:
$$\mathbf{A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T = U_k \Sigma_k V_k^T}$$
де $U_k \in \mathbb{R}^{m \times k}, \Sigma_k \in \mathbb{R}^{k \times k}, V_k \in \mathbb{R}^{n \times k}$.

* **Похибка наближення (Error bounds):**
  * Спектральна норма похибки: $\|A - A_k\|_2 = \sigma_{k+1}$ (дорівнює першому відкинутому сингулярному числу).
  * Норма Фробеніуса похибки: $\|A - A_k\|_F = \sqrt{\sum_{i=k+1}^r \sigma_i^2}$.

---

## 4. Практичні застосування SVD (Applications of SVD)

```text
                            Практичні застосування SVD
                                        |
         +--------------------+---------+---------+--------------------+
         |                    |                   |                    |
   Стиснення даних       Псевдообернення     Зниження шуму      Рекомендаційні системи
  (Image Compression)     Мура-Пенроуза       та NLP (LSA)       (Collaborative Filtering)
   A ≈ U_k Σ_k V_k^T      A^+ = V Σ^+ U^T     Видалення малих     Факторизація матриць
                                              сингулярних чисел   користувач-фільм (Netflix)
```

### 4.1. Стиснення зображень (Image Compression)
Чорно-біле зображення високої роздільності $m \times n$ займає $m \cdot n$ чисел.
* Збереження наближення рангу $k$: зберігаємо лише $U_k, \Sigma_k, V_k$, що вимагає:
  $$\text{Пам'ять} = k(m + n + 1) \ll m \cdot n$$
* Наприклад, для зображення $1000 \times 1000$ при $k = 50$: замість $1\,000\,000$ чисел зберігаємо всього $50 \cdot (1000 + 1000 + 1) \approx 100\,000$ чисел (стиснення в 10 разів із майже непомітною для людського ока втратою деталей!).

### 4.2. Псевдообернена матриця Мура-Пенроуза (Pseudoinverse $A^+$)
Через SVD псевдообернена матриця для **довільної прямокутної матриці** знаходиться миттєво:
$$\mathbf{A^+ = V \Sigma^+ U^T}$$
де $\Sigma^+ \in \mathbb{R}^{n \times m}$ утворюється заміною кожного ненульового елемента $\sigma_i$ на $1/\sigma_i$ та транспонуванням матриці:
$$\Sigma^+ = \text{diag}\left( \frac{1}{\sigma_1}, \frac{1}{\sigma_2}, \dots, \frac{1}{\sigma_r}, 0, \dots, 0 \right)$$
* Для будь-якої системи $A\mathbf{x} = \mathbf{b}$ вектор $\mathbf{x}_{\text{opt}} = A^+ \mathbf{b}$ є **єдиним розв'язком найменших квадратів з мінімальною евклідовою нормою** $\|\mathbf{x}\|_2$.

### 4.3. Латентно-семантичний аналіз (LSA / LSI in NLP)
Матриця «термін — документ» (Term-Document Matrix $A$):
* Рядки — слова словника, стовпці — статті / документи.
* Усічений розклад $A_k = U_k \Sigma_k V_k^T$ виділяє $k$ прихованих «тем» (латентних концепцій), групуючи синоніми та визначаючи семантичну схожість текстів.

---

## 📊 Зведена таблиця SVD (Summary Table)

| Поняття | Формула | Властивості |
|---|---|---|
| Сингулярний розклад | $A = U \Sigma V^T$ | $U^T U = I_m, \quad V^T V = I_n, \quad \sigma_i \ge 0$ |
| Сингулярні числа | $\sigma_i = \sqrt{\lambda_i(A^T A)}$ | Сортуються за спаданням: $\sigma_1 \ge \sigma_2 \ge \dots > 0$ |
| Найкраще наближення рангу $k$ | $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ | Похибка $\|A - A_k\|_2 = \sigma_{k+1}$ (Теорема Екхарта-Янга) |
| Псевдообернена матриця | $A^+ = V \Sigma^+ U^T$ | Розв'язує МНК для будь-якої $A \in \mathbb{R}^{m \times n}$ |

---

[⬅️ 12. Спектральна теорія та унітарні матриці](12-symmetric-unitary-matrices-and-spectral-decomposition.md) | [🏠 Головний зміст](index.md) | [⚡ Cheat Sheet](cheat-sheet.md)
