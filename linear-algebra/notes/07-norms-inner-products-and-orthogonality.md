# 07. Норми, відстані, скалярний добуток та ортогональність (Norms, Inner Products & Orthogonality)

[⬅️ 06. Заміна базису та оператори](06-change-of-basis-and-linear-transformations.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 08. Ортонормовані базиси, проекції та метод найменших квадратів ➡️](08-orthonormal-bases-projections-and-least-squares.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Скалярний добуток (Inner Product / Dot Product $\langle \mathbf{u}, \mathbf{v} \rangle$):** Алгебраїчна операція над векторами, що задає геометрію простору.
2. **Векторна норма (Vector Norm $\|\mathbf{v}\|$):** Міра довжини вектора ($L_1, L_2, L_\infty, L_p$).
3. **Метрика / Відстань (Distance $d(\mathbf{u}, \mathbf{v})$):** Довжина вектора різниці $\|\mathbf{u} - \mathbf{v}\|$.
4. **Нерівність Коші-Буняковського-Шварца (Cauchy-Schwarz Inequality):** Фундаментальне обмеження на скалярний добуток.
5. **Ортогональність векторів та підпросторів (Orthogonality & Orthogonal Complements):** Перпендикулярність векторів і взаємне доповнення підпросторів.

---

## 1. Скалярний добуток (Inner Product)

### 1.1. Аксіоматичне означення евклідового простору
Нехай $V$ — векторний простір над $\mathbb{R}$. **Скалярний добуток (Inner Product)** — це функція $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$, яка для будь-яких $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ та $c \in \mathbb{R}$ задовольняє **4 аксіоми**:
1. **Симетричність (Symmetry):** $\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$.
2. **Лінійність за першим аргументом (Linearity):**
   * $\langle \mathbf{u} + \mathbf{w}, \mathbf{v} \rangle = \langle \mathbf{u}, \mathbf{v} \rangle + \langle \mathbf{w}, \mathbf{v} \rangle$
   * $\langle c\mathbf{u}, \mathbf{v} \rangle = c \langle \mathbf{u}, \mathbf{v} \rangle$
3. **Невід'ємна визначеність (Positive Definiteness):**
   * $\langle \mathbf{v}, \mathbf{v} \rangle \ge 0$ для всіх $\mathbf{v} \in V$.
   * $\langle \mathbf{v}, \mathbf{v} \rangle = 0 \iff \mathbf{v} = \mathbf{0}$.

*(Для комплексних просторів $V$ над $\mathbb{C}$ — ермітів скалярний добуток: ермітова симетрія $\langle \mathbf{u}, \mathbf{v} \rangle = \overline{\langle \mathbf{v}, \mathbf{u} \rangle}$ та спряжена лінійність).*

### 1.2. Стандартний скалярний добуток (Dot Product в $\mathbb{R}^n$)
Для векторів $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$:
$$\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T \mathbf{v} = \sum_{i=1}^n u_i v_i = u_1 v_1 + u_2 v_2 + \dots + u_n v_n$$

---

## 2. Норми та відстані (Norms and Metrics)

### 2.1. Аксіоми норми (Norm Axioms)
**Норма (Norm $\|\cdot\|$):** Функція $\|\cdot\|: V \to \mathbb{R}$, яка зіставляє кожному вектору його «довжину» та задовольняє:
1. **Невід'ємність та строгість:** $\|\mathbf{v}\| \ge 0$, причому $\|\mathbf{v}\| = 0 \iff \mathbf{v} = \mathbf{0}$.
2. **Однорідність (Absolute Scalability):** $\|c\mathbf{v}\| = |c| \cdot \|\mathbf{v}\| \quad \forall c \in \mathbb{R}$.
3. **Нерівність трикутника (Triangle Inequality):** $\|\mathbf{u} + \mathbf{v}\| \le \|\mathbf{u}\| + \|\mathbf{v}\|$.

### 2.2. Індукована евклідова норма ($L_2$ Norm)
Будь-який скалярний добуток породжує природну норму:
$$\|\mathbf{v}\|_2 = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle} = \sqrt{\mathbf{v}^T \mathbf{v}} = \sqrt{\sum_{i=1}^n v_i^2}$$

### 2.3. Популярні $L_p$ норми в Data Science & Machine Learning

```text
    Одиничні сфери (Unit spheres in R^2: ||x|| <= 1):
           L1 (Манхеттен / Ромб)      L2 (Евклід / Коло)      L_inf (Чебишов / Квадрат)
                 ^ y                        ^ y                        ^ y
                 |   /\                     |   .-.                    |  +---+
                 |  /  \                    |  (   )                   |  |   |
               --+--\--/--> x             --+---`-'--> x             --+--+---+--> x
                 |   \/                     |                          |  +---+
```

1. **Манхеттенська норма ($L_1$ norm / Taxicab norm / Lasso regularization):**
   $$\|\mathbf{v}\|_1 = \sum_{i=1}^n |v_i| = |v_1| + |v_2| + \dots + |v_n|$$
   *(Сприяє розрідженості векторів / sparsity).*
2. **Евклідова норма ($L_2$ norm / Ridge regularization):**
   $$\|\mathbf{v}\|_2 = \sqrt{\sum_{i=1}^n v_i^2}$$
3. **Максимум-норма ($L_\infty$ norm / Chebyshev norm):**
   $$\|\mathbf{v}\|_\infty = \max_{1 \le i \le n} |v_i|$$
4. **Загальна $L_p$ норма ($p \ge 1$):**
   $$\|\mathbf{v}\|_p = \left( \sum_{i=1}^n |v_i|^p \right)^{1/p}$$

### 2.4. Відстань (Distance / Metric)
Відстань між векторами $\mathbf{u}$ та $\mathbf{v}$ визначається як норма їхньої різниці:
$$d(\mathbf{u}, \mathbf{v}) = \|\mathbf{u} - \mathbf{v}\|_2 = \sqrt{\sum_{i=1}^n (u_i - v_i)^2}$$

---

## 3. Фундаментальні нерівності та співвідношення

### 3.1. Нерівність Коші-Буняковського-Шварца (Cauchy-Schwarz Inequality)
Для будь-яких векторів $\mathbf{u}, \mathbf{v}$ у просторі зі скалярним добутком:
$$|\langle \mathbf{u}, \mathbf{v} \rangle| \le \|\mathbf{u}\| \cdot \|\mathbf{v}\| \iff |\mathbf{u}^T \mathbf{v}| \le \|\mathbf{u}\|_2 \|\mathbf{v}\|_2$$
Рівність досягається тоді і тільки тоді, коли вектори $\mathbf{u}$ та $\mathbf{v}$ є **лінійно залежними (колінеарними)**: $\mathbf{u} = c\mathbf{v}$.

### 3.2. Кут між векторами (Angle between vectors)
Завдяки нерівності Коші-Шварца значення $\frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{u}\| \|\mathbf{v}\|} \in [-1, 1]$, що дозволяє строго визначити **кут $\theta$ між двома ненульовими векторами**:
$$\cos\theta = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{u}\| \|\mathbf{v}\|} = \frac{\mathbf{u}^T \mathbf{v}}{\|\mathbf{u}\|_2 \|\mathbf{v}\|_2} \implies \theta = \arccos\left(\frac{\mathbf{u}^T \mathbf{v}}{\|\mathbf{u}\| \|\mathbf{v}\|}\right)$$
*(Це лежить в основі **косинусної подібності (Cosine Similarity)** в NLP та векторних базах даних).*

### 3.3. Теорема Піфагора та Закон паралелограма
1. **Теорема Піфагора (Pythagorean Theorem):**
   Якщо $\mathbf{u} \perp \mathbf{v}$ ($\langle \mathbf{u}, \mathbf{v} \rangle = 0$), то:
   $$\|\mathbf{u} + \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$$
2. **Тотожність паралелограма (Parallelogram Law):**
   $$\|\mathbf{u} + \mathbf{v}\|^2 + \|\mathbf{u} - \mathbf{v}\|^2 = 2\|\mathbf{u}\|^2 + 2\|\mathbf{v}\|^2$$
   *(Сума квадратів діагоналей паралелограма дорівнює сумі квадратів усіх чотирьох його сторін).*

---

## 4. Ортогональність векторів та підпросторів (Orthogonality)

### 4.1. Ортогональні вектори
Два вектори $\mathbf{u}, \mathbf{v} \in V$ називаються **ортогональними ($\mathbf{u} \perp \mathbf{v}$)**, якщо їхній скалярний добуток дорівнює нулю:
$$\mathbf{u} \perp \mathbf{v} \iff \langle \mathbf{u}, \mathbf{v} \rangle = 0 \iff \mathbf{u}^T \mathbf{v} = 0$$

* **Властивість:** Будь-яка множина попарно ортогональних ненульових векторів $\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$ є **лінійно незалежною**!

### 4.2. Ортогональні підпростори (Orthogonal Subspaces)
Два підпростори $V_1, V_2 \le V$ називаються **ортогональними ($V_1 \perp V_2$)**, якщо **кожен** вектор з $V_1$ перпендикулярний до **кожного** вектора з $V_2$:
$$\forall \mathbf{v}_1 \in V_1, \forall \mathbf{v}_2 \in V_2 \implies \langle \mathbf{v}_1, \mathbf{v}_2 \rangle = 0$$

### 4.3. Ортогональне доповнення (Orthogonal Complement $W^\perp$)
Нехай $W$ — підпростір у $V$. **Ортогональним доповненням $W^\perp$** називається множина всіх векторів у $V$, які перпендикулярні до кожного вектора з $W$:
$$W^\perp = \{\mathbf{v} \in V \mid \langle \mathbf{v}, \mathbf{w} \rangle = 0 \quad \forall \mathbf{w} \in W\}$$

* **Фундаментальні властивості:**
  1. $W^\perp \le V$ (є підпростором).
  2. $W \cap W^\perp = \{\mathbf{0}\}$.
  3. $(W^\perp)^\perp = W$.
  4. $\dim(W) + \dim(W^\perp) = \dim(V)$.
  5. Пряма сума: $V = W \oplus W^\perp$ (будь-який вектор $\mathbf{x} \in V$ однозначно розкладається на $\mathbf{x} = \mathbf{w} + \mathbf{w}^\perp$).

---

## 5. ⚠️ Типові помилки
* 🚨 **Плутанина між ортогональними підпросторами та перетином площин:** У 3D дві площини (наприклад, $xy$ та $yz$) перетинаються під кутом $90^\circ$, але вони **НЕ є ортогональними підпросторами**, оскільки їхній перетин — це пряма (вісь $y$), вздовж якої вектори належать обом площинам і не є перпендикулярними самі до себе! Для двох площин у 3D $\dim(W_1) + \dim(W_2) = 2+2=4 > 3$, що неможливо для ортогонального доповнення.

---

## 📊 Зведена таблиця (Summary Table)

| Поняття (Concept) | Математичний вираз | Геометричний зміст |
|---|---|---|
| Скалярний добуток | $\mathbf{u}^T \mathbf{v} = \|\mathbf{u}\| \|\mathbf{v}\| \cos\theta$ | Проекція векторів + кут |
| Евклідова норма | $\|\mathbf{v}\|_2 = \sqrt{\mathbf{v}^T \mathbf{v}}$ | Фізична довжина вектора |
| Нерівність Коші-Шварца | $|\mathbf{u}^T \mathbf{v}| \le \|\mathbf{u}\|_2 \|\mathbf{v}\|_2$ | Обмеження косинуса $[-1, 1]$ |
| Ортогональність | $\mathbf{u}^T \mathbf{v} = 0 \iff \theta = 90^\circ$ | Взаємна перпендикулярність |
| Ортогональне доповнення | $V = W \oplus W^\perp$ | Розклад простору на перпендикулярні частини |

---

[⬅️ 06. Заміна базису та оператори](06-change-of-basis-and-linear-transformations.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 08. Ортонормовані базиси, проекції та метод найменших квадратів ➡️](08-orthonormal-bases-projections-and-least-squares.md) | [⚡ Cheat Sheet](cheat-sheet.md)

