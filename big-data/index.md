# 🐘 Аналітика великих даних (Big Data Analytics with Databricks)

[🏠 Головна сторінка репозиторію](../README.md) | [⚡ Загальний Cheat Sheet курсу](cheat-sheet.md)

Цей розділ містить структуровані академічні конспекти, інженерні шпаргалки та банк екзаменаційних запитань (Q&A) з курсу **«Big Data Analytics Workflows with Databricks»** (UCU, 2026, викладач Андрій Музичук).

---

> **🔥 Швидкий доступ для підготовки до тестів:**
> Повна зведена інженерна шпаргалка з усіма термінами 3V, архітектурою Spark, типами трансформацій, форматами Avro/Parquet/Delta та патернами Data Quality доступна тут 👉 [**⚡ ULTIMATE BIG DATA CHEAT SHEET**](cheat-sheet.md)

---

## 🗺️ Навігація по лекціях курсу

| № | Назва лекції | Ключові теми та концепції (Core Concepts) | Шпаргалка | Q&A |
|---|---|---|---|---|
| ⭐ | **General Cheat Sheet** | **Зведена вижимка для тестів та співбесід по всьому курсу Big Data** | [**⚡ Cheat Sheet**](cheat-sheet.md) | — |
| **01** | **[Вступ до Big Data: Властивості, типи обробки та сховища даних](l1/index.md)** | Джерела даних, модель $3\text{V}$ (Volume, Velocity, Variety) + Veracity, Value, Variability; Batch vs Stream; NoSQL, Data Warehouse (Schema-on-Write), Data Lake (Schema-on-Read, Data Swamp), Lakehouse, платформа Databricks. | [⚡ Шпаргалка L1](l1/cheat-sheet.md) | [❓ Q&A L1](l1/questions.md) |
| **02** | **[Apache Spark: Архітектура, RDD, DataFrame та оптимізація](l2/index.md)** | Драйвер, воркери, екзекутори, слоти ($S = W \cdot C$), розрахунок хвиль задач; Application $\to$ Job $\to$ Stage $\to$ Task; RDD vs DataFrame; Narrow vs Wide трансформації; Shuffle проти Spill to Disk; керування пам'яттю (Execution vs Storage pool borrowing); Catalyst Optimizer, Tungsten, AQE, Broadcast Join. | [⚡ Шпаргалка L2](l2/cheat-sheet.md) | [❓ Q&A L2](l2/questions.md) |
| **03** | **[Формати файлів у Big Data та архітектура Delta Lake](l3/index.md)** | Рядковий (Row) проти колоночного (Columnar) збереження; Projection Pushdown, Predicate Pushdown (Data Skipping), колоночне стиснення; Avro (Kafka, JSON-схема), ORC (Hive), Parquet (Hadoop/Spark, Row Groups, Page); Delta Lake (`_delta_log`, ACID, Time Travel, Upsert/MERGE, OPTIMIZE, Z-ORDER); Iceberg та Hudi. | [⚡ Шпаргалка L3](l3/cheat-sheet.md) | [❓ Q&A L3](l3/questions.md) |
| **04** | **[Lakehouse: Моделювання даних, Medallion Architecture та Data Quality](l4/index.md)** | Bill Inmon (Top-down, 3NF EDW) проти Ralph Kimball (Bottom-up, Dimensional Modeling); нормалізація від UNF до 3NF та крок назад Кімболла; таблиці фактів (адитивні, напівадитивні, неадитивні), Grain; виміри Conformed, Role-playing, Junk; Схема «Зірка» vs «Сніжинка»; SCD Type 2 (Surrogate Key, інтервал $[effective, expiry)$, чому замало одного MERGE); Medallion Architecture (Bronze $\to$ Silver $\to$ Gold); Проблема понеділка: чому 0.07% браку ламає фінансовий `SUM()`, 4 стратегії (Ignore, Fail, Drop, Quarantine), рівняння закриття балансу ($\text{Received} = \text{Delivered} + \text{Rejected} + \text{Deduplicated}$) та Re-drive. | [⚡ Шпаргалка L4](l4/cheat-sheet.md) | [❓ Q&A L4](l4/questions.md) |

---

[🏠 Головна сторінка репозиторію](../README.md)
