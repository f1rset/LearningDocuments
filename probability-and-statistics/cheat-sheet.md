# ⚡ Інженерна шпаргалка: Теорія ймовірностей та математична статистика (Probability & Statistics Cheat Sheet)

[🏠 Головний зміст](index.md) | [📚 Повний курс лекцій](notes/01-probability-space-combinatorics-and-axioms.md)

---

## 📌 TL;DR — Квінтесенція курсу (Core Highlights)

* **Формула Байєса:** $P(H_k \mid A) = \frac{P(H_k)P(A \mid H_k)}{\sum P(H_i)P(A \mid H_i)}$. (Posterior $\propto$ Likelihood $\times$ Prior).
* **Лінійність сподівання:** $E[aX + bY] = aE[X] + bE[Y]$ (виконується **завжди**).
* **Дисперсія:** $\operatorname{Var}(X) = E[X^2] - (E[X])^2$, $\operatorname{Var}(aX+b) = a^2\operatorname{Var}(X)$.
* **Коваріація та кореляція:** $\operatorname{Cov}(X, Y) = E[XY] - \mu_X \mu_Y$, $\rho = \frac{\operatorname{Cov}}{\sigma_X \sigma_Y} \in [-1, 1]$. Незалежність $\implies \rho = 0$ (зворотне хибне, крім Гауссових!).
* **ЦГТ (Центральна гранична теорема):** $\frac{\bar{X}_n - \mu}{\sigma / \sqrt{n}} \xrightarrow{d} \mathcal{N}(0, 1)$ при $n \ge 30$.
* **Оцінка дисперсії:** $S^2 = \frac{1}{n-1}\sum (X_i - \bar{X})^2$ (поправка Бесселя $n-1$).
* **MLE (Метод максимальної правдоподібності):** $\hat{\theta} = \arg\max_\theta \sum \ln f(x_i; \theta)$.
* **Довірчий інтервал для середнього:** $\bar{X} \pm t_{n-1, 1-\alpha/2}\frac{S}{\sqrt{n}}$ (при невідомому $\sigma$).
* **Критерій $p$-value:** Якщо $p \le \alpha \implies$ відхилити $H_0$.
* **OLS Лінійна регресія:** $\hat{\boldsymbol{\beta}} = (X^T X)^{-1} X^T \mathbf{y}$, $\hat{\mathbf{y}} = H\mathbf{y}$.

---

## 1. Зведена таблиця розподілів ймовірностей (Probability Distributions Master Table)

### 1.1. Дискретні розподіли (Discrete Distributions)

| Розподіл | Позначення | PMF $P(X = k)$ | $E[X]$ | $\operatorname{Var}(X)$ | Опис / Застосування |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Бернуллі** | $\text{Bern}(p)$ | $p^k (1-p)^{1-k}, \; k \in \{0, 1\}$ | $p$ | $p(1-p)$ | 1 бінарна спроба |
| **Біноміальний** | $\text{Bin}(n, p)$ | $\binom{n}{k} p^k (1-p)^{n-k}, \; 0 \le k \le n$ | $np$ | $np(1-p)$ | Успіхи у вибірці з поверненням |
| **Геометричний** | $\text{Geom}(p)$ | $(1-p)^{k-1} p, \; k \in \{1, 2, \dots\}$ | $\frac{1}{p}$ | $\frac{1-p}{p^2}$ | Спроби до першого успіху (безпам'ятний) |
| **Від'ємний біноміальний**| $\text{NegBin}(r, p)$ | $\binom{k-1}{r-1} p^r (1-p)^{k-r}, \; k \ge r$ | $\frac{r}{p}$ | $\frac{r(1-p)}{p^2}$ | Спроби до $r$-го успіху |
| **Гіпергеометричний** | $\text{Hyper}(N, K, n)$| $\frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}$ | $n\frac{K}{N}$ | $n\frac{K}{N}(1-\frac{K}{N})\frac{N-n}{N-1}$ | Вибірка без повернення |
| **Пуассона** | $\text{Pois}(\lambda)$ | $\frac{\lambda^k e^{-\lambda}}{k!}, \; k \in \{0, 1, \dots\}$ | $\lambda$ | $\lambda$ | Рідкісні події, потік запитів |

### 1.2. Неперервні розподіли (Continuous Distributions)

| Розподіл | Позначення | Щільність PDF $f(x)$ | $E[X]$ | $\operatorname{Var}(X)$ | Опис / Застосування |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Рівномірний** | $U(a, b)$ | $\frac{1}{b-a}, \; x \in [a, b]$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | Рівноймовірні значення на відрізку |
| **Експоненційний** | $\text{Exp}(\lambda)$ | $\lambda e^{-\lambda x}, \; x \ge 0$ | $\frac{1}{\lambda}$ | $\frac{1}{\lambda^2}$ | Час очікування події (безпам'ятний) |
| **Нормальний** | $\mathcal{N}(\mu, \sigma^2)$| $\frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ | $\mu$ | $\sigma^2$ | Фундаментальний розподіл (дзвін) |
| **Стандартний нормальний** | $\mathcal{N}(0, 1)$ | $\frac{1}{\sqrt{2\pi}} e^{-z^2/2}$ | $0$ | $1$ | $Z = \frac{X - \mu}{\sigma}$ |
| **Логнормальний** | $\text{LogN}(\mu, \sigma^2)$ | $\frac{1}{x\sigma\sqrt{2\pi}} e^{-\frac{(\ln x - \mu)^2}{2\sigma^2}}$ | $e^{\mu + \sigma^2/2}$ | $e^{2\mu+\sigma^2}(e^{\sigma^2}-1)$ | Доходи, ціни акцій ($X > 0$) |
| **Гамма** | $\text{Gamma}(\alpha, \beta)$ | $\frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} e^{-\beta x}$ | $\frac{\alpha}{\beta}$ | $\frac{\alpha}{\beta^2}$ | Сума експоненційних величин |
| **Бета** | $\text{Beta}(\alpha, \beta)$ | $\frac{1}{B(\alpha, \beta)} x^{\alpha-1}(1-x)^{\beta-1}$ | $\frac{\alpha}{\alpha+\beta}$ | $\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$ | Розподіл імовірностей $x \in [0, 1]$ |

---

## 2. Вибіркові розподіли та зв'язки (Sampling Distributions Matrix)

| Статистика | Формула | Теоретичний розподіл | Призначення |
| :--- | :--- | :---: | :--- |
| **Нормальна $Z$** | $Z = \frac{\bar{X} - \mu}{\sigma / \sqrt{n}}$ | $\mathcal{N}(0, 1)$ | Середнє при відомому $\sigma$ |
| **Хі-квадрат $\chi^2$**| $\frac{(n-1)S^2}{\sigma^2} = \sum_{i=1}^n \frac{(X_i - \bar{X})^2}{\sigma^2}$ | $\chi^2(n - 1)$ | Дисперсія нормальної вибірки |
| **Стьюдента $t$** | $T = \frac{\bar{X} - \mu}{S / \sqrt{n}}$ | $t(n - 1)$ | Середнє при невідомому $\sigma$ |
| **Фішера $F$** | $F = \frac{S_1^2 / \sigma_1^2}{S_2^2 / \sigma_2^2}$ | $F(n_1 - 1, \; n_2 - 1)$ | Порівняння двох дисперсій (ANOVA) |

---

## 3. Зведена таблиця довірчих інтервалів (Confidence Intervals Lookup Table)

| Параметр | Умови | Формула $(1 - \alpha)$-довірчого інтервалу |
| :--- | :--- | :--- |
| **Середнє $\mu$** | $\sigma$ відоме | $\bar{X} \pm z_{1 - \alpha/2} \frac{\sigma}{\sqrt{n}}$ |
| **Середнє $\mu$** | $\sigma$ невідоме | $\bar{X} \pm t_{n-1, \; 1 - \alpha/2} \frac{S}{\sqrt{n}}$ |
| **Дисперсія $\sigma^2$** | Нормальна вибірка | $\left[ \frac{(n-1)S^2}{\chi^2_{n-1, \; 1-\alpha/2}}, \quad \frac{(n-1)S^2}{\chi^2_{n-1, \; \alpha/2}} \right]$ |
| **Різниця середніх $\mu_1 - \mu_2$** | Незалежні, $\sigma_1 = \sigma_2$ | $(\bar{X}_1 - \bar{X}_2) \pm t_{n_1+n_2-2, \; 1-\alpha/2} S_p \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}$ |
| **Різниця середніх (Парні)** | Залежні спостереження | $\bar{D} \pm t_{n-1, \; 1-\alpha/2} \frac{S_D}{\sqrt{n}}$ |
| **Частка успіхів $p$** | $n\hat{p} \ge 10$ (Вальд) | $\hat{p} \pm z_{1 - \alpha/2} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}}$ |

---

## 4. Гід по вибору статистичних критеріїв (Hypothesis Testing Guide)

```text
                  ЯКИЙ ТЕСТ ОБРАТИ?
                         |
       +-----------------+-----------------+
       |                                   |
   Числові дані                      Категорійні дані
       |                                   |
  Кількість груп?                    Таблиця спряженості?
  /    |    \                              |
 1     2     >2                            +---> Chi-Square (χ²) Test
 |     |      |                                  (Goodness-of-fit / Independence)
 |     |      +---> 1-Way ANOVA (Параметричний)
 |     |            Kruskal-Wallis (Непараметричний)
 |     |
 |     +---> Незалежні? 
 |           ├── ТАК: 2-Sample t-test / Welch (Нормальні) -> Mann-Whitney U (Ненорм.)
 |           └── НІ:  Paired t-test (Нормальні) -> Wilcoxon Signed-Rank (Ненорм.)
 |
 +---> 1-Sample t-test (при невідомому σ) / 1-Sample Z-test (при відомому σ)
```

---

## 5. Формули лінійної регресії (Linear Regression OLS Matrix Master Table)

| Поняття | Матрична формула | Скалярний випадок (1D) |
| :--- | :--- | :--- |
| **Оцінка коефіцієнтів $\hat{\boldsymbol{\beta}}$** | $\hat{\boldsymbol{\beta}} = (X^T X)^{-1} X^T \mathbf{y}$ | $\hat{\beta}_1 = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sum (x_i - \bar{x})^2} = r\frac{S_y}{S_x}$ |
| **Матриця-капелюх (Hat Matrix $H$)**| $H = X(X^T X)^{-1}X^T$ | $h_{ii} = \frac{1}{n} + \frac{(x_i - \bar{x})^2}{\sum (x_k - \bar{x})^2}$ |
| **Прогнозовані значення $\hat{\mathbf{y}}$** | $\hat{\mathbf{y}} = X\hat{\boldsymbol{\beta}} = H\mathbf{y}$ | $\hat{y}_i = \hat{\beta}_0 + \hat{\beta}_1 x_i$ |
| **Залишки моделі $\mathbf{e}$** | $\mathbf{e} = (I - H)\mathbf{y}$ | $e_i = y_i - \hat{y}_i$ |
| **Оцінка дисперсії помилок $s^2$** | $s^2 = \frac{\mathbf{e}^T \mathbf{e}}{n - p - 1} = \frac{\text{RSS}}{n - p - 1}$ | $s^2 = \frac{\sum e_i^2}{n - 2}$ |
| **Коваріація оцінок $\operatorname{Var}(\hat{\boldsymbol{\beta}})$** | $\operatorname{Var}(\hat{\boldsymbol{\beta}}) = s^2 (X^T X)^{-1}$ | $\operatorname{SE}(\hat{\beta}_1) = \frac{s}{\sqrt{\sum (x_i - \bar{x})^2}}$ |
| **Коефіцієнт детермінації $R^2$** | $R^2 = 1 - \frac{\text{RSS}}{\text{TSS}}$ | $R^2 = r_{XY}^2$ |
| **$F$-Статистика моделі** | $F = \frac{\text{ESS}/p}{\text{RSS}/(n - p - 1)}$ | $F = t_{\hat{\beta}_1}^2$ |

---

## 6. Типові підводні камені та інженерні правила (Rules of Thumb)

1. **Правило $n \ge 30$ для ЦГТ:** Для помірно симетричних розподілів вибіркове середнє стає нормальним вже при $n \ge 30$. Для сильно скошених розподілів може знадобитися $n \ge 100\text{--}500$.
2. **Кореляція $\neq$ Причинність:** Висока кореляція вказує на наявність зв'язку, але для доведення причинно-наслідкового впливу необхідний контрольований рандомізований A/B експеримент.
3. **Не потрапляйте у пастку $p$-Hacking:** Якщо провести 20 незалежних тестів при рівні $\alpha = 0.05$, хоча б один дасть «значущий» результат випадково з імовірністю $64\%$. Використовуйте поправку Бонферроні ($\alpha / m$) або Benjamini-Hochberg (FDR).

---

[🏠 Головний зміст](index.md) | [📚 Повний курс лекцій](notes/01-probability-space-combinatorics-and-axioms.md)

