# 🎲 Теорія ймовірностей та математична статистика: Повний курс-конспект (Probability & Statistics Complete Course)

[🏠 Головна сторінка репозиторію](../../README.md) | [⚡ Ultimate Cheat Sheet](cheat-sheet.md)

Цей курс містить вичерпні, глибокі та академічно строгі матеріали з **Теорії ймовірностей, математичної статистики, теорії оцінювання, перевірки статистичних гіпотез та регресійного аналізу**. Всі терміни та назви теорем повністю дублюються **англійською мовою**.

---

> **🔥 Швидкий доступ:** Повна інженерна шпаргалка з усіма таблицями розподілів, вибіркових статистик, довірчих інтервалів, матрицею вибору критеріїв та формулами регресії доступна тут 👉 [**⚡ ULTIMATE PROBABILITY & STATISTICS CHEAT SHEET**](cheat-sheet.md)

---

## 🗺️ Навігація по модулях курсу (Modules & Topics)

| № | Тема / Модуль | Ключові концепції (Core Concepts) | Англійські терміни | Посилання |
|---|---|---|---|---|
| ⭐ | **Cheat Sheet** | **Таблиця розподілів, вибіркові розподіли, довірчі інтервали, критерії, регресія OLS** | **Ultimate Prob & Stat Cheat Sheet** | [**⚡ Відкрити Cheat Sheet**](cheat-sheet.md) |
| **01** | **Простір подій, комбінаторика та аксіоматика** | Алгебра подій, аксіоми Колмогорова ($\Omega, \mathcal{F}, P$), класична/геометрична ймовірність, комбінаторика ($P_n, A_n^k, C_n^k$), включення-виключення. | *Probability Space, Kolmogorov Axioms, Combinatorics, Permutations, Combinations, Inclusion-Exclusion* | [Читати модуль 01 ➡️](notes/01-probability-space-combinatorics-and-axioms.md) |
| **02** | **Умовна ймовірність, Байєс та незалежність** | Умовна ймовірність $P(A\|B)$, формула повної ймовірності, формула Байєса (апріорне/апостеріорне), незалежність подій, схема Бернуллі. | *Conditional Probability, Law of Total Probability, Bayes' Theorem, Independence, Bernoulli Trials* | [Читати модуль 02 ➡️](notes/02-conditional-probability-bayes-and-independence.md) |
| **03** | **Дискретні випадкові величини та розподіли** | PMF, CDF, розподіли: Бернуллі, Біноміальний, Геометричний (безпам'ятність), Від'ємний біноміальний, Гіпергеометричний, Пуассона. | *Discrete Random Variables, PMF, CDF, Bernoulli, Binomial, Geometric, Poisson Distribution* | [Читати модуль 03 ➡️](notes/03-discrete-random-variables-and-distributions.md) |
| **04** | **Неперервні випадкові величини та щільності** | PDF $f(x)$, CDF $F(x)$, Рівномірний, Експоненційний (безпам'ятний, hazard rate), Нормальний $\mathcal{N}(\mu, \sigma^2)$, $Z$-шкала, правило $3\sigma$, Гамма, Бета, Логнормальний, заміна змінних. | *Continuous Random Variables, PDF, Uniform, Exponential, Normal (Gaussian), $3\sigma$ Rule, Gamma, Beta* | [Читати модуль 04 ➡️](notes/04-continuous-random-variables-and-densities.md) |
| **05** | **Числові характеристики та моменти** | Математичне сподівання $E[X]$ (лінійність, LOTUS), дисперсія $\operatorname{Var}(X)$, стандартне відхилення $\sigma$, асиметрія (Skewness), ексцес (Kurtosis), MGF $M_X(t)$, характеристична функція $\varphi_X(t)$. | *Expected Value, Linearity, Variance, Higher Moments, Skewness, Kurtosis, MGF, Characteristic Function* | [Читати модуль 05 ➡️](notes/05-expectation-variance-and-moments.md) |
| **06** | **Багатовимірні розподіли та коваріація** | Спільна щільність $f_{X,Y}$, маргінали, умовне сподівання $E[Y\|X]$, коваріація, кореляція Пірсона $\rho$, коваріаційна матриця $\Sigma \succeq 0$, багатовимірний Гауссів вектор. | *Joint & Marginal Distributions, Law of Total Expectation, Covariance, Correlation, Covariance Matrix, Multivariate Normal* | [Читати модуль 06 ➡️](notes/06-multivariate-distributions-covariance-and-joint-pdf.md) |
| **07** | **Граничні теореми: ЗВЧ та ЦГТ** | Нерівності Маркова, Чебишова, Єнсена; режими збіжності ($\text{a.s.}, P, L_2, d$), теорема Слуцького, Закон великих чисел (WLLN, SLLN), Центральна гранична теорема (ЦГТ), Муавр-Лаплас. | *Probability Inequalities (Markov, Chebyshev, Jensen), Modes of Convergence, Law of Large Numbers, Central Limit Theorem* | [Читати модуль 07 ➡️](notes/07-limit-theorems-law-of-large-numbers-and-clt.md) |
| **08** | **Описова статистика та вибіркові розподіли** | Генеральна сукупність vs Вибірка, вибіркове середнє $\bar{X}$, виправлена дисперсія $S^2$ ($n-1$), ECDF, теорема Гливенка-Кантеллі, розподіли $\chi^2(k), t(k), F(d_1, d_2)$, теорема Фішера. | *Sample Statistics, Bessel's Correction, ECDF, Chi-Squared, Student's $t$, Fisher's $F$, Fisher's Theorem* | [Читати модуль 08 ➡️](notes/08-descriptive-statistics-and-sampling-distributions.md) |
| **09** | **Точкове оцінювання, метод моментів та MLE** | Незміщеність, спроможність, ефективність, MSE (Bias-Variance tradeoff), інформація Фішера $I(\theta)$, нерівність Рао-Крамера (CRLB), метод моментів (MoM), максимальна правдоподібність (MLE), MAP. | *Point Estimation, Unbiasedness, Consistency, CRLB, Method of Moments, Maximum Likelihood (MLE), MAP* | [Читати модуль 09 ➡️](notes/09-point-estimation-mle-and-method-of-moments.md) |
| **10** | **Довірчі інтервали та інтервальне оцінювання** | Рівень довіри $1-\alpha$, $Z$-інтервал, $t$-інтервал Стьюдента, $\chi^2$-інтервал дисперсії, різниця середніх двох вибірок (пулова $S_p^2$, Велч, парний), інтервали Вальда/Вільсона для частоти $p$, Бутстреп. | *Confidence Intervals, Pivot Method, $t$-Interval, Difference of Means, Proportion CIs (Wilson), Bootstrap* | [Читати модуль 10 ➡️](notes/10-confidence-intervals.md) |
| **11** | **Перевірка статистичних гіпотез** | $H_0$ vs $H_1$, помилки I та II роду ($\alpha, \beta$), потужність тесту $1-\beta$, $p$-value, $Z$-тест, $t$-тести, $F$-тест, 1-Way ANOVA, $\chi^2$-критерій згоди та незалежності, тест Манна-Вітні, множинне тестування (Бонферроні, FDR). | *Hypothesis Testing, Type I & II Errors, Statistical Power, $p$-Value, $t$-Test, ANOVA, Chi-Square Test, Mann-Whitney, FDR* | [Читати модуль 11 ➡️](notes/11-hypothesis-testing.md) |
| **12** | **Лінійна регресія та кореляційний аналіз** | Проста та множинна регресія, матричний МНК $\hat{\boldsymbol{\beta}} = (X^T X)^{-1}X^T\mathbf{y}$, матриця-капелюх $H$, теорема Гаусса-Маркова (BLUE), $R^2$ та $R^2_{\text{adj}}$, $t$-/$F$-тести коефіцієнтів, діагностика залишків ($\text{VIF}$, гомоскедастичність). | *Linear Regression, Ordinary Least Squares (OLS), Hat Matrix, Gauss-Markov Theorem, $R^2$, Residual Diagnostics, Multicollinearity* | [Читати модуль 12 ➡️](notes/12-linear-regression-and-correlation-analysis.md) |

---

[🏠 Головна сторінка репозиторію](../../README.md)

