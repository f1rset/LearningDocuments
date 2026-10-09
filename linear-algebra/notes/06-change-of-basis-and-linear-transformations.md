# 06. Заміна базису, матриці переходу та лінійні оператори (Change of Basis & Linear Transformations)

[⬅️ 05. Базиси та 4 підпростори](05-bases-dimension-rank-and-four-subspaces.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Норми, відстані та ортогональність ➡️](07-norms-inner-products-and-orthogonality.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Матриця переходу (Change of Basis Matrix / Transition Matrix $P$):** Перерахунок координат векторів між різними базисами.
2. **Лінійне відображення / оператор (Linear Transformation / Linear Map $T$):** Відображення, що зберігає операції додавання та множення на скаляр.
3. **Матриця лінійного відображення (Matrix Representation of $T$):** Стовпці матриці — це образи базисних векторів.
4. **Матриця оператора в різних базисах (Similar Matrices):** Подібність матриць $B = P^{-1} A P$.

---

## 1. Заміна базису та матриця переходу (Change of Basis)

### 1.1. Координатний вектор
Нехай у просторі $V$ задано базис $\mathcal{B} = \{\mathbf{v}_1, \mathbf{v}_2, \dots, \mathbf{v}_n\}$.
Будь-який вектор $\mathbf{x} \in V$ записується як:
$$\mathbf{x} = x_1 \mathbf{v}_1 + x_2 \mathbf{v}_2 + \dots + x_n \mathbf{v}_n \iff [\mathbf{x}]_{\mathcal{B}} = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}$$

### 1.2. Матриця переходу (Transition Matrix)
Нехай є два базиси в $V$:
* Старий базис (Old basis): $\mathcal{B} = \{\mathbf{v}_1, \dots, \mathbf{v}_n\}$
* Новий базис (New basis): $\mathcal{B}' = \{\mathbf{u}_1, \dots, \mathbf{u}_n\}$

Виразимо вектори нового базису через вектори старого базису:
$$\mathbf{u}_j = \sum_{i=1}^n p_{ij} \mathbf{v}_i \iff [\mathbf{u}_j]_{\mathcal{B}} = \begin{bmatrix} p_{1j} \\ p_{2j} \\ \vdots \\ p_{nj} \end{bmatrix}$$

**Матриця переходу від $\mathcal{B}$ до $\mathcal{B}'$ (Transition Matrix $P_{\mathcal{B} \to \mathcal{B}'}$):**
Це матриця, стовпцями якої є координати нових базисних векторів у старому базисі:
$$P = P_{\mathcal{B} \to \mathcal{B}'} = \begin{bmatrix} [\mathbf{u}_1]_{\mathcal{B}} & [\mathbf{u}_2]_{\mathcal{B}} & \dots & [\mathbf{u}_n]_{\mathcal{B}} \end{bmatrix}$$

### 1.3. Формула перетворення координат вектора
$$[\mathbf{x}]_{\mathcal{B}} = P_{\mathcal{B} \to \mathcal{B}'} [\mathbf{x}]_{\mathcal{B}'} \iff \mathbf{[\mathbf{x}]_{\mathcal{B}'} = P_{\mathcal{B} \to \mathcal{B}'}^{-1} [\mathbf{x}]_{\mathcal{B}}}$$

> **⚠️ Зверніть увагу на напрямок:**
> * Матриця $P$ складається з координат **нових** векторів у **старому** базисі.
> * Проте для знаходження **нових** координат вектора $[\mathbf{x}]_{\mathcal{B}'}$ множать на **$P^{-1}$**!

---

## 2. Лінійні відображення та оператори (Linear Transformations)

### 2.1. Означення лінійного відображення
Відображення $T: V \to W$ між векторними просторами над полем $\mathbb{F}$ називається **лінійним (linear transformation / linear map)**, якщо для всіх $\mathbf{u}, \mathbf{v} \in V$ та $c \in \mathbb{F}$ виконуються дві властивості:
1. **Адитивність (Additivity):** $T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v})$.
2. **Однорідність (Homogeneity):** $T(c\mathbf{u}) = c T(\mathbf{u})$.
*(Разом: $T(c\mathbf{u} + d\mathbf{v}) = c T(\mathbf{u}) + d T(\mathbf{v})$).*

Якщо $V = W$, відображення $T: V \to V$ називається **лінійним оператором (linear operator / endomorphism)**.

### 2.2. Властивості лінійних відображень:
* $T(\mathbf{0}_V) = \mathbf{0}_W$.
* $T(-\mathbf{v}) = -T(\mathbf{v})$.
* **Ядро лінійного відображення (Kernel):** $\ker(T) = \{\mathbf{v} \in V \mid T(\mathbf{v}) = \mathbf{0}\} \le V$.
* **Образ лінійного відображення (Image / Range):** $\text{im}(T) = \{T(\mathbf{v}) \mid \mathbf{v} \in V\} \le W$.
* **Теорема про ранг і дефект для операторів:** $\dim(\ker(T)) + \dim(\text{im}(T)) = \dim(V)$.

---

## 3. Матриця лінійного відображення (Matrix of a Linear Transformation)

### 3.1. Побудова матриці оператора
Нехай $T: V \to W$, де:
* Базис простору $V$: $\mathcal{B}_V = \{\mathbf{v}_1, \dots, \mathbf{v}_n\}$, $\dim(V) = n$.
* Базис простору $W$: $\mathcal{B}_W = \{\mathbf{w}_1, \dots, \mathbf{w}_m\}$, $\dim(W) = m$.

Щоб знайти матрицю $[T]_{\mathcal{B}_V \to \mathcal{B}_W} \in \mathbb{R}^{m \times n}$:
1. Застосовуємо оператор $T$ до кожного базисного вектора вхідного простору: $T(\mathbf{v}_1), \dots, T(\mathbf{v}_n)$.
2. Розкладаємо кожен отриманий образ за базисом $\mathcal{B}_W$: $[T(\mathbf{v}_j)]_{\mathcal{B}_W}$.
3. Записуємо координатні стовпці у матрицю:
   $$[T]_{\mathcal{B}_V \to \mathcal{B}_W} = \begin{bmatrix} [T(\mathbf{v}_1)]_{\mathcal{B}_W} & [T(\mathbf{v}_2)]_{\mathcal{B}_W} & \dots & [T(\mathbf{v}_n)]_{\mathcal{B}_W} \end{bmatrix}$$

Тоді обчислення образу будь-якого вектора зводиться до матричного множення:
$$[T(\mathbf{x})]_{\mathcal{B}_W} = [T]_{\mathcal{B}_V \to \mathcal{B}_W} \cdot [\mathbf{x}]_{\mathcal{B}_V}$$

---

## 4. Матриця оператора при зміні базису (Similar Matrices)

### 4.1. Формула зміни базису оператора
Нехай $T: V \to V$ — лінійний оператор, і задано два базиси в $V$: $\mathcal{B}$ та $\mathcal{B}'$.
Нехай $P = P_{\mathcal{B} \to \mathcal{B}'}$ — матриця переходу від $\mathcal{B}$ до $\mathcal{B}'$.
Нехай $A = [T]_{\mathcal{B}}$ — матриця оператора у старому базисі $\mathcal{B}$.
Нехай $B = [T]_{\mathcal{B}'}$ — матриця оператора у новому базисі $\mathcal{B}'$.

```text
                  [x]_B'   --- B --->   [T(x)]_B'
                    |                      ^
                    | P                    | P^-1
                    v                      |
                  [x]_B    --- A --->   [T(x)]_B
```

Формула зв'язку між матрицями в різних базисах:
$$\mathbf{B = P^{-1} A P}$$

### 4.2. Подібні матриці (Similar Matrices)
Дві квадратні матриці $A$ та $B$ називаються **подібними ($A \sim B$)**, якщо існує невироджена матриця $P$ така, що:
$$B = P^{-1} A P$$

#### Інваріанти подібних матриць (не залежать від вибору базису):
1. **Визначник:** $\det(B) = \det(P^{-1} A P) = \det(P^{-1})\det(A)\det(P) = \det(A)$.
2. **Слід (Trace):** $\text{tr}(B) = \text{tr}(A)$.
3. **Ранг:** $\text{rank}(B) = \text{rank}(A)$.
4. **Характеристичний многочлен:** $p_B(\lambda) = \det(B - \lambda I) = \det(A - \lambda I) = p_A(\lambda)$.
5. **Власні значення (Eigenvalues):** Множина власних значень (спектр) матриць $A$ та $B$ є ідентичною!

---

## 5. ⚠️ Типові помилки (Common Pitfalls)
* 🚨 **Плутанина між $P$ та $P^{-1}$ у формулі подібності:** Запам'ятайте схему: щоб застосувати старий оператор $A$ до нових координат $[\mathbf{x}]_{\mathcal{B}'}$, ми спочатку переходимо в старий базис ($P [\mathbf{x}]_{\mathcal{B}'}$), діємо оператором $A$ ($A P [\mathbf{x}]_{\mathcal{B}'}$), а потім повертаємо результат у новий базис ($P^{-1} A P$).
* 🚨 **Нелінійні перетворення:** Зсув на вектор ($T(\mathbf{x}) = \mathbf{x} + \mathbf{x}_0$ при $\mathbf{x}_0 \neq \mathbf{0}$) **НЕ є лінійним відображенням** (бо $T(\mathbf{0}) \neq \mathbf{0}$). Такі відображення є афінними (affine transformations).

---

## 📊 Зведена таблиця (Summary Table)

| Концепція | Формула | Опис |
|---|---|---|
| Матриця переходу | $P = [[\mathbf{u}_1]_{\mathcal{B}} \dots [\mathbf{u}_n]_{\mathcal{B}}]$ | Стовпці — нові вектори в старому базисі |
| Перерахунок вектора | $[\mathbf{x}]_{\mathcal{B}'} = P^{-1} [\mathbf{x}]_{\mathcal{B}}$ | Зміна координат вектора |
| Матриця відображення | $A = [[T(\mathbf{v}_1)]_{\mathcal{W}} \dots [T(\mathbf{v}_n)]_{\mathcal{W}}]$ | Стовпці — образи базисних векторів |
| Зміна базису оператора | $B = P^{-1} A P$ | $A$ та $B$ — подібні матриці (Similar) |

---

[⬅️ 05. Базиси та 4 підпростори](05-bases-dimension-rank-and-four-subspaces.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 07. Норми, відстані та ортогональність ➡️](07-norms-inner-products-and-orthogonality.md) | [⚡ Cheat Sheet](cheat-sheet.md)
