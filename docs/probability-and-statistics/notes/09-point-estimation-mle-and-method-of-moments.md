# 09. Точкове оцінювання параметрів та метод MLE (Point Estimation, Method of Moments & MLE)

[⬅️ Попередня: 08. Описова статистика та вибіркові розподіли](08-descriptive-statistics-and-sampling-distributions.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Довірчі інтервали ➡️](10-confidence-intervals.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Точкова оцінка ($\hat{\theta}$):** Статистика для оцінки невідомого параметра генеральної сукупності $\theta$.
2. **Критерії якості оцінок:** Незміщеність ($E[\hat{\theta}] = \theta$), Спроможність ($\hat{\theta} \xrightarrow{P} \theta$), Ефективність (MSE = $\operatorname{Var} + \text{Bias}^2$).
3. **Інформація Фішера ($I(\theta)$) та нерівність Рао-Крамера (CRLB):** Теоретична нижня межа дисперсії будь-якої незміщеної оцінки.
4. **Метод моментів (Method of Moments - MoM):** Прирівнювання вибіркових моментів до теоретичних.
5. **Метод максимальної правдоподібності (Maximum Likelihood Estimation - MLE):** Максимізація логарифмічної правдоподібності $\ell(\theta) = \ln L(\theta)$ та асимптотична ефективність.
6. **Байєсівське оцінювання (Bayesian Inference) та MAP:** Спряжені апріорні розподіли та зв'язок регуляризації з апріорним знанням.

---

## 1. Критерії якості точкових оцінок (Properties of Estimators)

Нехай $X_1, \dots, X_n \sim f(x; \theta)$ — вибірка, а $\hat{\theta} = T(X_1, \dots, X_n)$ — точкова оцінка параметра $\theta$.

### 1.1. Незміщеність (Unbiasedness)
Оцінка $\hat{\theta}$ називається **незміщеною**, якщо її математичне сподівання точно дорівнює істинному значенню параметра:
$$E[\hat{\theta}] = \theta \iff \operatorname{Bias}(\hat{\theta}) = E[\hat{\theta}] - \theta = 0$$
* *Приклад:* $\bar{X}$ є незміщеною оцінкою $\mu$, а $S^2 = \frac{1}{n-1}\sum(X_i-\bar{X})^2$ — незміщеною оцінкою $\sigma^2$.

### 1.2. Спроможність (Consistency)
Оцінка $\hat{\theta}_n$ називається **спроможною**, якщо при збільшенні розміру вибірки вона збігається до $\theta$ за ймовірністю:
$$\hat{\theta}_n \xrightarrow{P} \theta \iff \lim_{n \to \infty} P(|\hat{\theta}_n - \theta| < \varepsilon) = 1, \quad \forall \varepsilon > 0$$

### 1.3. Середньоквадратична похибка (Mean Squared Error - MSE)
$$\operatorname{MSE}(\hat{\theta}) = E\left[ (\hat{\theta} - \theta)^2 \right] = \operatorname{Var}(\hat{\theta}) + \left(\operatorname{Bias}(\hat{\theta})\right)^2$$
*(Компроміс зміщення та дисперсії / Bias-Variance Tradeoff: іноді злегка зміщена оцінка має значно меншу дисперсію і дає менший MSE)*.

```text
    Низький Bias, Висока Var        Високий Bias, Низька Var       Ідеально (Низький Bias, Низька Var)
           .  *  .                          . . * .                        . * .
         *   θ   *                        *  .  . *                       *  θ  *
           .  *  .                          . . * .                        . * .
                                               θ
```

---

## 2. Інформація Фішера та нерівність Рао-Крамера (CRLB)

### 2.1. Інформація Фішера (Fisher Information $I(\theta)$)
Кількість інформації про параметр $\theta$, яку несе одне випадкове спостереження:
$$I(\theta) = E\left[ \left( \frac{\partial \ln f(X; \theta)}{\partial \theta} \right)^2 \right] = - E\left[ \frac{\partial^2 \ln f(X; \theta)}{\partial \theta^2} \right]$$
* Для вибірки розміру $n$ сумарна інформація дорівнює $I_n(\theta) = n I(\theta)$.

### 2.2. Нерівність Рао-Крамера (Cramér-Rao Lower Bound - CRLB)
> **👑 Теорема Рао-Крамера:**
> Дисперсія **будь-якої** незміщеної оцінки $\hat{\theta}$ не може бути меншою за величину, обернену до повної інформації Фішера:
> $$\operatorname{Var}(\hat{\theta}) \ge \frac{1}{n I(\theta)}$$

* Якщо $\operatorname{Var}(\hat{\theta}) = \frac{1}{n I(\theta)}$, оцінка $\hat{\theta}$ називається **ефективною (Efficient)** або UMVUE (Uniformly Minimum-Variance Unbiased Estimator).

---

## 3. Метод моментів (Method of Moments - MoM)

### 3.1. Принцип методу
Прирівнюємо теоретичні моменти розподілу $\alpha_k(\theta) = E[X^k]$ до відповідних вибіркових емпіричних моментів $m_k = \frac{1}{n}\sum_{i=1}^n X_i^k$:

$$\begin{cases}
E[X] = \frac{1}{n}\sum_{i=1}^n X_i = \bar{X} \\
E[X^2] = \frac{1}{n}\sum_{i=1}^n X_i^2 \\
\dots
\end{cases}$$

* **Плюси:** Простий у розрахунках, завжди дає спроможні оцінки.
* **Мінуси:** Оцінки часто є неефективними (мають більшу дисперсію, ніж MLE).

---

## 4. Метод максимальної правдоподібності (Maximum Likelihood Estimation - MLE)

### 4.1. Функція правдоподібності (Likelihood Function)
Для фіксованої вибірки даних $\mathbf{x} = (x_1, \dots, x_n)$ функція правдоподібності розглядається як функція від параметра $\theta$:
$$L(\theta; \mathbf{x}) = \prod_{i=1}^n f(x_i; \theta)$$

### 4.2. Логарифмічна функція правдоподібності (Log-Likelihood)
Оскільки логарифм є строго монотонною функцією, точка максимуму $L(\theta)$ збігається з точкою максимуму $\ell(\theta) = \ln L(\theta)$:
$$\ell(\theta) = \ln L(\theta; \mathbf{x}) = \sum_{i=1}^n \ln f(x_i; \theta)$$

* **Рівняння правдоподібності (Score Equation):**
  $$\frac{\partial \ell(\theta)}{\partial \theta} = 0 \iff \sum_{i=1}^n \frac{\partial \ln f(x_i; \theta)}{\partial \theta} = 0$$

### 4.3. Асимптотичні властивості MLE
При $n \to \infty$ оцінка максимальної правдоподібності $\hat{\theta}_{\text{MLE}}$ має ідеальні статистичні властивості:
1. **Спроможність:** $\hat{\theta}_{\text{MLE}} \xrightarrow{P} \theta$.
2. **Асимптотична нормальність:**
   $$\sqrt{n}(\hat{\theta}_{\text{MLE}} - \theta) \xrightarrow{d} \mathcal{N}\left( 0, \; \frac{1}{I(\theta)} \right)$$
3. **Асимптотична ефективність:** Досягає нижньої межі Рао-Крамера!
4. **👑 Інваріантність (Invariance Property):** Якщо $\hat{\theta}$ — оцінка MLE для $\theta$, то для будь-якої функції $g(\theta)$ оцінкою MLE є $g(\hat{\theta})$.

---

## 5. Байєсівське оцінювання та MAP (Bayesian Estimation & MAP)

У Байєсівському підході невідомий параметр $\theta$ вважається не фіксованим числом, а **випадковою величиною** з відомим апріорним розподілом $p(\theta)$.

* **Теорема Байєса для параметрів:**
  $$p(\theta \mid \mathbf{x}) = \frac{p(\mathbf{x} \mid \theta) p(\theta)}{p(\mathbf{x})} = \frac{L(\theta; \mathbf{x}) p(\theta)}{\int L(\theta; \mathbf{x}) p(\theta)\,d\theta} \propto L(\theta; \mathbf{x}) \cdot p(\theta)$$
  $$\text{Posterior} \propto \text{Likelihood} \times \text{Prior}$$

### 5.1. Точкова оцінка максимуму апостеріорної ймовірності (MAP)
$$\hat{\theta}_{\text{MAP}} = \arg\max_\theta p(\theta \mid \mathbf{x}) = \arg\max_\theta \left[ \ell(\theta) + \ln p(\theta) \right]$$
* *Зв'язок із Machine Learning:* 
  * Нормальний апріорний розподіл $p(\theta) \sim \mathcal{N}(0, \tau^2) \implies L_2$-регуляризація (Ridge).
  * Розподіл Лапласа $p(\theta) \propto e^{-\lambda|\theta|} \implies L_1$-регуляризація (Lasso).

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Зміщеність вибіркової дисперсії в MLE для Гаусса:**
   Оцінка MLE дисперсії для нормального розподілу $\hat{\sigma}^2_{\text{MLE}} = \frac{1}{n}\sum(x_i - \bar{x})^2$ є **зміщеною** (має знаменник $n$, а не $n-1$). Хоча асимптотично при $n \to \infty$ зміщення зникає.
2. **Максимум на межі області визначення (Uniform):**
   Для вибірки $X_i \sim U(0, \theta)$ похідна правдоподібності не перетворюється на нуль. Максимум досягається на межі: $\hat{\theta}_{\text{MLE}} = \max(X_1, \dots, X_n) = X_{(n)}$.
3. **Плутанина між MLE та MAP:** MLE шукає максимум правдоподібності даних $L(\theta)$, тоді як MAP враховує апріорну інформацію $p(\theta)$, зсуваючи оцінку в бік апріорного середнього.

---

[⬅️ Попередня: 08. Описова статистика та вибіркові розподіли](08-descriptive-statistics-and-sampling-distributions.md) | [🏠 Головний зміст](index.md) | [Наступна тема: 10. Довірчі інтервали ➡️](10-confidence-intervals.md) | [⚡ Cheat Sheet](cheat-sheet.md)
