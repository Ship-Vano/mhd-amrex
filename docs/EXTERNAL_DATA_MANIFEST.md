# Manifest внешнего архива `Elements`

Состояние чтения: 2026-08-21. Том смонтирован как `/Volumes/Elements`.
Во время аудита файлы на внешнем диске не изменялись и не копировались в Git.
Полный архив результатов имеет размер около 60 GB, поэтому в репозитории
хранятся только пути, метаданные и контрольные суммы ключевых объектов.

## 1. Найденные копии исходного кода

### `/Volumes/Elements/MHD2D`

- remote: `https://github.com/Ship-Vano/MHD2D.git`;
- Git HEAD: `3023cb0003038cbd8b38c77aaf24123d5815d106` от 2025-06-08;
- рабочее дерево не является чистым snapshot: на FAT/exFAT изменены modes и
  окончания строк, а `MHDSolver2D.cpp` содержит незакоммиченные изменения;
- после игнорирования modes/CRLF рабочий `MHDSolver2D.cpp` отличается от
  baseline `9d0f60e…` тремя строками, относящимися к осесимметричному
  MHD-blast (`t_end`, `Bz`, индекс EMF), а не к Brio–Wu;
- `debug/InputData` содержит постановку `taskType=3`, не Brio–Wu;
- `debug/OutputData` занимает около 207 MB и содержит 77 файлов.

Ключевые debug inputs:

| Файл | Размер | SHA-256 |
|---|---:|---|
| `debug/InputData/mesh.txt` | 932498 B | `cade070f32623c8b48088b882d15509b2fd59465c5d175613f72de52b7dec6c6` |
| `debug/InputData/netConfig.json` | 153 B | `bf4a62fa12b957ece6520e6243f204d78f1e7a72d8026046f7bf4e77e4cade33` |
| `debug/InputData/solverConfig.json` | 265 B | `d4ae8840faacb3b88506afc4c2132c2cb86f182bedf8dad205a5b6075810a482` |

Сетка имеет 9992 узла и 19592 треугольника. Конфиг задаёт `t_end=0.01`,
но рабочий исходник для `taskType=3` переопределяет его значением `0.003`.

### `/Volumes/Elements/MHD_program/MHD2D`

- Git HEAD: `db0a8043d6413b1743849499bdbb24df5f209ed1` от 2025-01-30;
- `InputData` и часть `OutputData` добавлены в index, но никогда не вошли в
  commit; остальной worktree также содержит изменения;
- `solverConfig.json` задаёт `taskType=2`, `t_end=0.9`, mesh `mesh3` — это
  Alfvén wave, а не Brio–Wu.

| Файл | Сетка/размер | SHA-256 |
|---|---:|---|
| `InputData/mesh.txt` | 2225 узлов, 3552 элемента | `3bcda581c0a0da6518cfc2bface99bb93d0a9f5ca51b695a8d29667c8c5d42ed` |
| `InputData/mesh3.txt` | 67553 узла, 134144 элемента | `7b4dccba32a65f067f99291c7dbf7fa9814a62950f562ee8f4c5072295c5689d` |
| `InputData/netConfig.json` | 156 B | `f5109ac34e337127ec7081a4fdec06755e0d962cf18cb8397827424a6122ae4a` |
| `InputData/solverConfig.json` | 184 B | `0068a9c20c88f67e393e8784b9d5a4cc9cfd3a34937b2c746b3273d4acaf00e1` |

## 2. Mesher

Каталог `/Volumes/Elements/MHD_program/meshGen` является отдельным Git
репозиторием:

- remote: `https://github.com/Ship-Vano/meshGen`;
- HEAD: `e30e2ffc16806948fb8efe869469abc51ebe8925` от 2025-01-30;
- фактически используется Python-библиотека **Netgen**, а не Gmsh;
- tracked: `main.py`, `structured.py`;
- untracked: `eqlit.py`, `mesh.txt`, `equilateral_mesh.txt`;
- рабочий `main.py` изменён относительно commit: для Brio–Wu `maxh` заменён
  с `0.02` на `0.009`, число refinements — с `5` на `2`.

| Файл | SHA-256 |
|---|---|
| рабочий `main.py` | `b40c0503fcd17bc4a668b791bc1dd52a05fab28c412ef6e710a223c94b409140` |
| `main.py` из commit | `88dbef4dad2f15e2a133a7454f21e03e238f9457539fd0c24c5d45815b70357a` |
| `structured.py` | `2044f6ba669db6ef24d198225569d58678612055402c00e43f3699d9284661cd` |
| `eqlit.py` | `c1a7a2a110e211692518bfa8024148fd102f2b50d5f67abae27c9a95d6d1c200` |
| созданный `mesh.txt` | `3bcda581c0a0da6518cfc2bface99bb93d0a9f5ca51b695a8d29667c8c5d42ed` |

Последний hash совпадает с `MHD_program/MHD2D/InputData/mesh.txt`, что даёт
надёжную связь mesher → input. Связи с `MHDresults/BrioWu/mesh5.bin` пока нет.

## 3. Архив результатов

Корень: `/Volumes/Elements/MHDresults`, общий размер около 60 GB. Среди
каталогов присутствуют Brio–Wu, две Alfvén wave серии, Orszag–Tang разных
размеров, rotor/rotating cylinder, magnetic loop, Cartesian/cylindrical Sod,
MHD blast и `diplomTests`.

### Brio–Wu

Каталог `/Volumes/Elements/MHDresults/BrioWu`:

- размер около 23 GB;
- `mesh5.vtu`: 32623 точки и 64576 треугольных ячеек;
- 1548 файлов `tmpres_*.vtu`;
- основной ряд: шаги `0,10,...,15460`;
- `tmpres_151.vtu` создан раньше основного ряда и является финалом отдельного
  короткого запуска; его нельзя включать в длинную временную серию;
- VTU содержит одну девятикомпонентную cell array `elemUs`;
- physical time в VTU не записан;
- рядом нет JSON-конфига, console log или commit manifest.

| Файл | SHA-256 |
|---|---|
| `mesh5.bin` | `da40448583fc1f50ef46314fbb82434f306c2d5cd78118425554225b5ff197d2` |
| `mesh5.vtu` | `6431e4c2493a051da40170d5d936f4d825b3d1bbf142439e0e3099427397ee6d` |
| `tmpres_0.vtu` | `8245e15a7b6ba094d92d20af012a03c640612929efdbbc3c4240c3bd0b0d69e0` |
| `tmpres_15460.vtu` | `44c24ed316d84893dc37c44f4a3b99decb296dc13ed9ed3c82196b853184fa58` |
| `brioWuDensity.0000.png` | `668c2c30151760eaed5a36b4220d72e5774cefd6581fd9607ec1ec2141a86b79` |

Основной ряд создан 2025-02-16. Последний Git commit `MHD2D` до этой даты —
`233f7de…` от 2025-02-14; в нём Brio–Wu использует CFL `0.4`, `gamma=2`,
`Bx=0.75`, free-flow boundaries и default `t_end=0.1`. Это только временная
корреляция, а не доказательство ревизии расчёта: worktree и JSON не сохранены.

### Частично восстановленные журналы ВКР

В `MHDresults/diplomTests` найдены одиночные CPU console logs:

| Case | Время из console | SHA-256 console |
|---|---:|---|
| Alfvén N=32 | 19.89687109 s | `a1b8e89702e1b19dc922adbfa4a323498a849eecbbdb9bac117b51f365b394be` |
| Alfvén N=64 | 164.4318438 s | `6940b2484df8e7e32f9b513b0b57919fbed24592d80fe7204bdb33ca115eb5d4` |
| Alfvén N=128 | 765.0536875 s | `f03420642d4782d4ec719aef3e37f422411990ebe55143fd950f6396dd464dfd` |
| magnetic loop N=100 | 1919.081 s | `3aa0648362bbee01a3fbf1c9e0931e212c50115496ada96f5d33e5e4dcdc7ca6` |

Alfvén N=128 согласуется с округлённым CPU-временем 765.1 s в ВКР. Это
частичная первичная запись, но не полноценный benchmark: один повтор, нет
machine manifest/flags/commit и GPU-log. Соответствующие 538.7, 342.13 и
215.39 s для CPU+GPU в найденных console logs отсутствуют.

## 4. Правила дальнейшего использования

1. Не добавлять 60 GB результатов в обычный Git.
2. Перед использованием копировать только выбранные inputs/logs и сохранять
   source path, размер, timestamp и SHA-256.
3. Не считать имя каталога доказательством case/config/commit.
4. Для Brio–Wu сначала написать VTU inspector, проверить порядок девяти
   компонентов `elemUs`, geometry, initial state и физическое время.
5. Исторический ряд использовать как `legacy_unattributed`, пока конфиг и
   точная ревизия не будут восстановлены либо независимо подтверждены.
6. Для нового benchmark выполнить расчёт заново; журналы ВКР использовать
   только как историческое сопоставление.
