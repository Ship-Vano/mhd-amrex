# Состав Git-репозитория

Документ фиксирует политику initial commit. Точный состав каждого коммита всё
равно проверяется командами `git status --short` и `git diff --cached`.

## Версионируемые файлы

- `.gitignore`, `AGENTS.md`, `CMakeLists.txt`, `README.md`;
- `src/**` — исходный код решателя и численных ядер;
- `inputs/*.json` — воспроизводимые постановки задач;
- `scripts/*.py` — анализ и визуализация;
- `tests/*.cpp` — исходники проверок;
- `docs/*.md`, `docs/*.tex` — требования, аудит и исходник отчёта;
- `docs/figures/*.png` — небольшие канонические рисунки, необходимые отчёту.

## Локальные или воспроизводимые артефакты

В Git не добавляются:

- `build*/`, локальные исполняемые файлы и CMake cache;
- AMReX plotfile `plt*` и автоматически созданные `.old.*` копии;
- ADIOS2 BP, HDF5, VTK/VTU и каталог `results/`;
- `tests/verify`, тестовые CSV/PNG;
- скомпилированный `docs/report.pdf` и LaTeX auxiliary files;
- профили Nsight, coverage, Python/editor/OS cache.

Сырые benchmark-данные могут быть слишком велики для Git. Если они хранятся
внешне, в репозиторий обязательно добавляются их manifest/checksum, скрипт
агрегации и итоговые компактные таблицы из `benchmarks/summary/`.

## Перед коммитом

```bash
git status --short --ignored
git diff --cached --check
git diff --cached --stat
```

Нельзя использовать `git add -f` для обхода правил без объяснения причины в
commit message или сопроводительном отчёте.
