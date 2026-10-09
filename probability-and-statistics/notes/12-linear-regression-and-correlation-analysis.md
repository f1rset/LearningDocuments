# 12. Лінійна регресія та кореляційний аналіз (Linear Regression & Correlation Analysis)

[⬅️ Попередня: 11. Перевірка статистичних гіпотез](11-hypothesis-testing.md) | [🏠 Головний зміст](index.md) | [⚡ Cheat Sheet](cheat-sheet.md)

---

## 🎯 Ключові концепції (Core Concepts)
1. **Лінійна модель регресії ($Y = X\boldsymbol{\beta} + \boldsymbol{\varepsilon}$):** Моделювання математичного сподівання відгуку залежно від предикторів.
2. **Метод найменших квадратів (OLS):** Матрична формула $\hat{\boldsymbol{\beta}} = (X^T X)^{-1} X^T \mathbf{y}$ та проєкційна матриця-капелюх $H$.
3. **Теорема Гаусса-Маркова (BLUE):** Оптимальність OLS серед усіх лінійних незміщених оцінок.
4. **Коефіцієнт детермінації ($R^2$ та $R^2_{\text{adj}}$):** Частка дисперсії, пояснена моделлю.
5. **Статистичний висновок для коефіцієнтів:** $t$-тести для окремих предикторів та $F$-тест загальної значущості регресії.
6. **Діагностика залишків (Residual Diagnostics):** Гетероскедастичність, мультиколінеарність ($\text{VIF}$) та нормальність залишків (Q-Q plot).

---

## 1. Постановка задачі лінійної регресії (Linear Model)

### 1.1. Проста лінійна регресія (Simple Linear Regression - SLR)
$$y_i = \beta_0 + \beta_1 x_i + \varepsilon_i, \quad i = 1, \dots, n$$
* Оцінки МНК для 1D випадку:
  $$\hat{\beta}_1 = \frac{\operatorname{Cov}(X, Y)}{\operatorname{Var}(X)} = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sum (x_i - \bar{x})^2} = r_{XY} \frac{S_Y}{S_X}, \qquad \hat{\beta}_0 = \bar{y} - \hat{\beta}_1 \bar{x}$$

### 1.2. Множинна лінійна регресія в матричній формі (Multiple Linear Regression - MLR)
$$\mathbf{y} = X \boldsymbol{\beta} + \boldsymbol{\varepsilon}$$
де:
$$\mathbf{y} = \begin{pmatrix} y_1 \\ \vdots \\ y_n \end{pmatrix} \in \mathbb{R}^n, \quad 
X = \begin{pmatrix} 1 & x_{11} & \dots & x_{1p} \\ \vdots & \vdots & \ddots & \vdots \\ 1 & x_{n1} & \dots & x_{np} \end{pmatrix} \in \mathbb{R}^{n \times (p+1)}, \quad 
\boldsymbol{\beta} = \begin{pmatrix} \beta_0 \\ \beta_1 \\ \vdots \\ \beta_p \end{pmatrix} \in \mathbb{R}^{p+1}, \quad 
\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_1 \\ \vdots \\ \varepsilon_n \end{pmatrix} \in \mathbb{R}^n$$

---

## 2. Метод найменших квадратів (Ordinary Least Squares - OLS)

### 2.1. Виведення оцінки $\hat{\boldsymbol{\beta}}$
Мінімізуємо суму квадратів залишків (Residual Sum of Squares - RSS):
$$\operatorname{RSS}(\boldsymbol{\beta}) = \|\mathbf{y} - X\boldsymbol{\beta}\|^2 = (\mathbf{y} - X\boldsymbol{\beta})^T (\mathbf{y} - X\boldsymbol{\beta}) \to \min_{\boldsymbol{\beta}}$$
Диференціюємо за вектором $\boldsymbol{\beta}$ та прирівнюємо до нуля (Нормальні рівняння / Normal Equations):
$$\nabla_{\boldsymbol{\beta}} \operatorname{RSS} = -2 X^T (\mathbf{y} - X\boldsymbol{\beta}) = \mathbf{0} \iff X^T X \boldsymbol{\beta} = X^T \mathbf{y}$$

> **👑 Формула OLS-оцінки:**
> $$\hat{\boldsymbol{\beta}} = (X^T X)^{-1} X^T \mathbf{y}$$

* **Прогнозовані значення (Fitted Values $\hat{\mathbf{y}}$) та матриця капелюх ($H$ / Hat Matrix):**
  $$\hat{\mathbf{y}} = X \hat{\boldsymbol{\beta}} = \underbrace{X (X^T X)^{-1} X^T}_{H} \mathbf{y} = H \mathbf{y}$$
* **Залишки (Residuals $\mathbf{e}$):**
  $$\mathbf{e} = \mathbf{y} - \hat{\mathbf{y}} = (I - H) \mathbf{y}$$
  *(Залишки ортогональні до простору стовпців матриці $X$: $X^T \mathbf{e} = \mathbf{0}$)*.

```text
               y (Вектор спостережень)
              /|
             / |
            /  | e = y - y_hat (Залишки: перпендикуляр)
           /   |
          *----+-------------> Простір стовпців Col(X)
         0     y_hat = H y (Ортогональна проекція)
```

---

## 3. Теорема Гаусса-Маркова та класичні припущення (Gauss-Markov)

### 3.1. Умови класичної лінійної моделі
1. **Лінійність:** Модель лінійна за невідомими параметрами $\boldsymbol{\beta}$.
2. **Строга екзогенність:** $E[\boldsymbol{\varepsilon} \mid X] = \mathbf{0}$ (предиктори не корелюють із шумом).
3. **Гомоскедастичність та відсутність автокореляції:**
   $$\operatorname{Var}(\boldsymbol{\varepsilon} \mid X) = \sigma^2 I_n \iff \operatorname{Var}(\varepsilon_i) = \sigma^2, \quad \operatorname{Cov}(\varepsilon_i, \varepsilon_j) = 0 \; (i \neq j)$$
4. **Повний стовпчиковий ранг:** $\operatorname{rank}(X) = p + 1 \le n$ (відсутність мультиколінеарності).
5. *(Для статистичних тестів)* **Нормальність шуму:** $\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \sigma^2 I_n)$.

> **👑 Теорема Гаусса-Маркова:**
> При виконанні умов 1–4 оцінка МНК $\hat{\boldsymbol{\beta}}$ є **BLUE (Best Linear Unbiased Estimator)** — найкращою лінійною незміщеною оцінкою (вона має найменшу дисперсію серед усіх лінійних незміщених оцінок).

---

## 4. Якість підгонки моделі: Розклад дисперсії та $R^2$

$$\underbrace{\sum_{i=1}^n (y_i - \bar{y})^2}_{\text{TSS (Total)}} = \underbrace{\sum_{i=1}^n (\hat{y}_i - \bar{y})^2}_{\text{ESS (Explained)}} + \underbrace{\sum_{i=1}^n e_i^2}_{\text{RSS (Residual)}}$$

* **Коефіцієнт детермінації ($R^2$):**
  $$R^2 = \frac{\text{ESS}}{\text{TSS}} = 1 - \frac{\text{RSS}}{\text{TSS}} \in [0, 1]$$
  *(Частка варіації залежної змінної $Y$, що пояснюється предикторами $X$)*.
* **Скоригований коефіцієнт детермінації (Adjusted $R^2$):**
  $$R^2_{\text{adj}} = 1 - \frac{\text{RSS} / (n - p - 1)}{\text{TSS} / (n - 1)} = 1 - (1 - R^2)\frac{n - 1}{n - p - 1}$$
  *(Штрафує за додавання непотрібних змінних до моделі)*.

---

## 5. Статистичний висновок у регресії (Inference in Regression)

* **Незміщена оцінка дисперсії помилок ($\hat{\sigma}^2$ / $s^2$):**
  $$s^2 = \frac{\text{RSS}}{n - p - 1} = \frac{\sum e_i^2}{n - p - 1}$$
* **Коваріаційна матриця оцінок коефіцієнтів:**
  $$\operatorname{Var}(\hat{\boldsymbol{\beta}}) = \sigma^2 (X^T X)^{-1} \implies \widehat{\operatorname{Var}}(\hat{\boldsymbol{\beta}}) = s^2 (X^T X)^{-1}$$

### 5.1. $t$-Тест для окремого коефіцієнта (Significance of Predictor)
$$H_0: \beta_j = 0 \quad \text{vs} \quad H_1: \beta_j \neq 0$$
$$T = \frac{\hat{\beta}_j - 0}{\text{SE}(\hat{\beta}_j)} \sim t(n - p - 1), \quad \text{де } \text{SE}(\hat{\beta}_j) = s \sqrt{[(X^T X)^{-1}]_{jj}}$$

### 5.2. $F$-Тест загальної значущості моделі (Overall Model $F$-Test)
$$H_0: \beta_1 = \beta_2 = \dots = \beta_p = 0 \quad \text{vs} \quad H_1: \exists \beta_j \neq 0$$
$$F = \frac{\text{ESS} / p}{\text{RSS} / (n - p - 1)} = \frac{R^2 / p}{(1 - R^2) / (n - p - 1)} \sim F(p, \; n - p - 1)$$

---

## 6. Діагностика залишків та порушення припущень (Diagnostics)

```text
    Ідеальні залишки (Гомоскедастичність)        Гетероскедастичність (Віяло / Нестала дисперсія)
      e ^  . * . . * .                             e ^        . * . *
        |  * . * . * .                               |      . * . *
        +------------------> y_hat                   +------------------> y_hat
        |  . * . . * .                               |      . * . *
                                                     |        . * . *
```

1. **Гетероскедастичність (Heteroscedasticity):** Дисперсія помилок $\operatorname{Var}(\varepsilon_i)$ залежить від $x$. 
   * *Наслідки:* Оцінки залишаються незміщеними, але стандартні помилки $\text{SE}(\hat{\beta})$ викривляються ($t$-тести хибні).
   * *Лікування:* Робастні стандартні помилки Вайта (HC1/HC3), зважений МНК (WLS), логарифмування $Y$.
2. **Мультиколінеарність (Multicollinearity):** Сильна лінійна залежність між самими предикторами.
   * *Діагностика:* Коефіцієнт здуття дисперсії $\text{VIF}_j = \frac{1}{1 - R_j^2} > 5\text{--}10$.
   * *Лікування:* Видалення зайвих змінних, регуляризація (Ridge/Lasso), PCA.
3. **Ненормальність залишків:** Перевіряється за допомогою Q-Q plot та тесту Жарка-Бера (Jarque-Bera).

---

## ⚠️ Підводні камені та типові помилки (Pitfalls)

1. **Кореляція НЕ означає причинно-наслідковий зв'язок (Correlation $\neq$ Causation):**
   Високий коефіцієнт $R^2$ свідчить про наявність статистичної асоціації, але не доводить, що зміна $X$ спричиняє зміну $Y$ (можливий спільний латентний фактор — Confounder).
2. **Парадокс додавання змінних до $R^2$:**
   Звичайний $R^2$ **завжди зростає** при додаванні будь-якої нової змінної, навіть випадкового шуму. Тому для порівняння моделей різної складності слід використовувати **$R^2_{\text{adj}}$**, AIC або BIC.
3. **Екстраполяція за межі діапазону даних:**
   Прогнозування за допомогою регресії далеко за межами діапазону вихідних значень $X$ призводить до катастрофічних помилок.

---

[⬅️ Попередня: 11. Перевірка статистичних гіпотез](11-hypothesis-testing.md) | [🏠 Головний зміст](index.md) | [⚡ Cheat Sheet](cheat-sheet.md)
