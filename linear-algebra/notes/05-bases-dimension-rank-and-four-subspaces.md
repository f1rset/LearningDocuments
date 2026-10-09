# 05. Базиси векторних просторів, ранг та 4 фундаментальні підпростори (Bases, Rank & The Four Fundamental Subspaces)

[⬅️ 04. Лінійні простори](04-vector-spaces-subspaces-and-spans.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 06. Заміна базису та лінійні оператори ➡️](06-change-of-basis-and-linear-transformations.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Базис (Basis):** Мінімальна породжуюча та максимальна лінійно незалежна множина векторів.
2. **Розмірність (Dimension / $\dim(V)$):** Кількість векторів у будь-якому базисі простору $V$.
3. **Координати та ізоморфізм (Coordinates & Isomorphism):** Будь-який скінченновимірний простір $\dim(V)=n$ ізоморфний $\mathbb{R}^n$.
4. **Чотири фундаментальні підпростори матриці (The Four Fundamental Subspaces):** Простір стовпців, простір рядків, ядро та ліве ядро.
5. **Основна теорема лінійної алгебри (The Fundamental Theorem of Linear Algebra) та теорема про ранг і дефект (Rank-Nullity Theorem).**

---

## 1. Базис та Розмірність (Basis and Dimension)

### 1.1. Означення базису
Множина векторів $\mathcal{B} = \{\mathbf{v}_1, \mathbf{v}_2, \dots, \mathbf{v}_n\} \subset V$ називається **базисом (basis)** векторного простору $V$, якщо:
1. Вектори $\{\mathbf{v}_1, \dots, \mathbf{v}_n\}$ **лінійно незалежні** (Linearly Independent).
2. Вектори $\{\mathbf{v}_1, \dots, \mathbf{v}_n\}$ **породжують простір $V$** ($\text{span}(\mathcal{B}) = V$).

> **Ключова властивість базису:** Кожен вектор $\mathbf{x} \in V$ виражається через базисні вектори **єдиним чином**:
> $$\mathbf{x} = c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \dots + c_n \mathbf{v}_n$$
> Числа $(c_1, c_2, \dots, c_n)$ називаються **координатами вектора $\mathbf{x}$ у базисі $\mathcal{B}$**:
> $$[\mathbf{x}]_{\mathcal{B}} = \begin{bmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{bmatrix} \in \mathbb{R}^n$$

### 1.2. Розмірність простору (Dimension)
* **Розмірність $\dim(V)$** — це кількість векторів у будь-якому базисі простору $V$.
* Якщо $\dim(V) = n$, то:
  * Будь-які $n$ лінійно незалежних векторів утворюють базис.
  * Будь-які $n$ векторів, що породжують $V$, утворюють базис.
  * Будь-яка множина з $> n$ векторів є лінійно залежною.
  * Будь-яка множина з $< n$ векторів не може породити $V$.

### 1.3. Ізоморфізм векторних просторів (Isomorphism)
Два векторні простори $V$ та $W$ над $\mathbb{F}$ називаються **ізоморфними ($V \cong W$)**, якщо існує бієктивне лінійне відображення $T: V \to W$.
* **Теорема про ізоморфізм:** Будь-який $n$-вимірний дійсний векторний простір $V$ ($\dim(V)=n$) є **ізоморфним до $\mathbb{R}^n$**:
  $$T: V \to \mathbb{R}^n, \quad \mathbf{v} \mapsto [\mathbf{v}]_{\mathcal{B}}$$
  *(Це дозволяє звести роботу з многочленами, матрицями чи абстрактними функціями до звичайних числових векторів в $\mathbb{R}^n$).*

---

## 2. Чотири фундаментальні підпростори матриці (The Four Fundamental Subspaces)

Для будь-якої матриці $A \in \mathbb{R}^{m \times n}$ рангу $\text{rank}(A) = r$ існують 4 природні підпростори:

```text
               Простір входів R^n                    Простір виходів R^m
        +-------------------------------+     +-------------------------------+
        |                               |     |                               |
        |   Простір рядків Row(A)       |     |   Простір стовпців Col(A)     |
        |   dim = r                     | ===>|   dim = r                     |
        |   (Базис: ненульові рядки REF)|     |   (Базис: pivot стовпці A)    |
        |                               |     |                               |
        +-------------------------------+     +-------------------------------+
        |  |_ (Ортогональне доповнення) |     |  |_ (Ортогональне доповнення) |
        +-------------------------------+     +-------------------------------+
        |                               |     |                               |
        |   Ядро Null(A)                |     |   Ліве ядро Null(A^T)         |
        |   dim = n - r                 |     |   dim = m - r                 |
        |   (Ax = 0)                    |     |   (A^T y = 0)                 |
        |                               |     |                               |
        +-------------------------------+     +-------------------------------+
```

### 2.1. Означення та розмірності

1. **Простір стовпців (Column Space / Image / Range / $\text{Col}(A)$ або $\text{im}(A)$):**
   * $\text{Col}(A) = \{\mathbf{y} \in \mathbb{R}^m \mid \mathbf{y} = A\mathbf{x} \text{ для деякого } \mathbf{x} \in \mathbb{R}^n\} \le \mathbb{R}^m$.
   * Лінійна оболонка стовпців: $\text{Col}(A) = \text{span}\{\mathbf{a}_1, \dots, \mathbf{a}_n\}$.
   * **Розмірність:** $\dim(\text{Col}(A)) = \text{rank}(A) = r$.
   * **Як знайти базис:** Стовпці вихідної матриці $A$, які відповідають **стовпцям з ведучими елементами (pivots)** в $\text{RREF}(A)$.

2. **Простір рядків (Row Space / $\text{Row}(A) = \text{Col}(A^T)$):**
   * $\text{Row}(A) = \{\mathbf{v} \in \mathbb{R}^n \mid \mathbf{v} = A^T \mathbf{y}\} \le \mathbb{R}^n$.
   * Лінійна оболонка рядків матриці $A$.
   * **Розмірність:** $\dim(\text{Row}(A)) = \text{rank}(A) = r$.
   * **Як знайти базис:** Ненульові рядки матриці $\text{REF}(A)$ або $\text{RREF}(A)$.

3. **Ядро / Простір нулів (Nullspace / Kernel / $\text{Null}(A)$ або $\ker(A)$):**
   * $\text{Null}(A) = \{\mathbf{x} \in \mathbb{R}^n \mid A\mathbf{x} = \mathbf{0}\} \le \mathbb{R}^n$.
   * **Розмірність (Дефект / Nullity):** $\dim(\text{Null}(A)) = n - r$.
   * **Як знайти базис:** Розв'язати $A\mathbf{x} = \mathbf{0}$ через RREF, виразивши головні змінні через $n-r$ вільних змінних. Вектори при вільних змінних утворюють фундаментальну систему розв'язків (ФСР).

4. **Ліве ядро (Left Nullspace / $\text{Null}(A^T)$):**
   * $\text{Null}(A^T) = \{\mathbf{y} \in \mathbb{R}^m \mid A^T \mathbf{y} = \mathbf{0} \iff \mathbf{y}^T A = \mathbf{0}^T\} \le \mathbb{R}^m$.
   * **Розмірність:** $\dim(\text{Null}(A^T)) = m - r$.
   * **Як знайти базис:** Розв'язати однорідну систему $A^T \mathbf{y} = \mathbf{0}$.

---

## 3. Фундаментальні теореми лінійної алгебри

### 3.1. Теорема про ранг і дефект (Rank-Nullity Theorem)
Для будь-якої матриці $A \in \mathbb{R}^{m \times n}$:
$$\mathbf{\text{rank}(A) + \dim(\text{Null}(A)) = n}$$
$$\text{Кількість pivot-стовпців} + \text{Кількість вільних змінних} = \text{Загальна кількість стовпців } n$$

### 3.2. Основна теорема лінійної алгебри: Ортогональність підпросторів
Підпростори матриці є **взаємно ортогональними доповненнями (orthogonal complements)** у відповідних просторах:
1. **В $\mathbb{R}^n$:**
   $$\text{Row}(A) \perp \text{Null}(A), \quad \text{Row}(A) \oplus \text{Null}(A) = \mathbb{R}^n \implies \text{Null}(A) = (\text{Row}(A))^\perp$$
2. **В $\mathbb{R}^m$:**
   $$\text{Col}(A) \perp \text{Null}(A^T), \quad \text{Col}(A) \oplus \text{Null}(A^T) = \mathbb{R}^m \implies \text{Null}(A^T) = (\text{Col}(A))^\perp$$

> **💡 Геометричний наслідок:**
> Будь-який вектор $\mathbf{x} \in \mathbb{R}^n$ розкладається на ортогональні складові: $\mathbf{x} = \mathbf{x}_{\text{row}} + \mathbf{x}_{\text{null}}$, де $\mathbf{x}_{\text{row}} \in \text{Row}(A)$ та $\mathbf{x}_{\text{null}} \in \text{Null}(A)$. При цьому $A\mathbf{x} = A\mathbf{x}_{\text{row}} + A\mathbf{x}_{\text{null}} = A\mathbf{x}_{\text{row}} + \mathbf{0} = A\mathbf{x}_{\text{row}}$. Матриця $A$ відображає підпростір $\text{Row}(A)$ **взаємно однозначно (бієктивно)** на $\text{Col}(A)$!

---

## 4. ⚠️ Типові помилки при пошуку базисів
* 🚨 **Базис $\text{Col}(A)$ беруть з RREF:** Базисні вектори простору стовпців треба брати **з початкової матриці $A$** (номери стовпців, де в RREF є pivots), а не з самої матриці RREF, оскільки рядкові операції Гаусса змінюють простір стовпців!
* 🚨 **Базис $\text{Row}(A)$ можна брати з RREF:** Рядкові операції НЕ змінюють простір рядків, тому ненульові рядки RREF є чудовим базисом для $\text{Row}(A)$.

---

## 📊 Зведена таблиця 4 підпросторів (Summary of 4 Subspaces)

| Підпростір (Subspace) | Позначення | Де живе? | Розмірність | Базис |
|---|---|---|---|---|
| **Простір стовпців (Column Space)** | $\text{Col}(A)$ | $\mathbb{R}^m$ | $r$ | Pivot-стовпці початкової матриці $A$ |
| **Ядро / Простір нулів (Nullspace)** | $\text{Null}(A)$ | $\mathbb{R}^n$ | $n - r$ | Вектори ФСР системи $A\mathbf{x}=\mathbf{0}$ |
| **Простір рядків (Row Space)** | $\text{Row}(A)$ | $\mathbb{R}^n$ | $r$ | Ненульові рядки $\text{RREF}(A)$ |
| **Ліве ядро (Left Nullspace)** | $\text{Null}(A^T)$ | $\mathbb{R}^m$ | $m - r$ | Вектори ФСР системи $A^T \mathbf{y}=\mathbf{0}$ |

---

[⬅️ 04. Лінійні простори](04-vector-spaces-subspaces-and-spans.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 06. Заміна базису та лінійні оператори ➡️](06-change-of-basis-and-linear-transformations.md) | [⚡ Cheat Sheet](cheat-sheet.md)
