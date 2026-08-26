# T00 baseline manifest

Дата фиксации: 2026-08-26. Этот manifest описывает вход T00 и не является
benchmark-артефактом.

## Репозиторий и рабочее дерево

| Поле | Значение |
|---|---|
| Репозиторий | `mhd-amrex` |
| Входной commit | `2e5089f7ba931d240d702d512c8f248d11f39882` |
| Состояние перед T00 | clean (`git status --short` пуст) |
| AMReX | 25.01, `FetchContent`, фиксированный `GIT_TAG` |
| Компилятор | Homebrew GCC 15.2.0 (`/opt/homebrew/opt/gcc/bin/g++-15`) |
| CMake | 4.3.3 |
| MPI | Open MPI 5.0.9 |
| OpenMP | включён, версия 4.5 |
| Платформа | macOS 15.7.3, 10 logical CPUs |

## Воспроизведённая команда

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DMHD_MPI=ON -DMHD_OPENMP=ON -DMHD_HDF5=OFF
cmake --build build -j 4
./build/mhd2d_verify briowu
./build/mhd2d_verify alfven 32
```

На входном commit команды завершились успешно. Наблюдаемые диагностические
выводы: для standalone Brio--Wu -- 488 шагов и `max|divB| = 0`; для
Alfvén `N=32` -- 80 шагов и историческая норма `max|divB| = 3.706e-13`.
Это smoke-проверки без CTest, сохранённого run manifest или независимого
эталона; они не верифицируют HLLD-ветви, AMR-консервацию, GPU-путь либо
производительность.

## Внешние источники и их статус

| Источник | Статус | Основание |
|---|---|---|
| legacy source | кандидат: `/Users/ivansamanov/Documents/MHD2D` @ `9d0f60e…` | [LEGACY_BASELINE.md](LEGACY_BASELINE.md) |
| внешние mesh/input/result | только инвентаризированы | [EXTERNAL_DATA_MANIFEST.md](EXTERNAL_DATA_MANIFEST.md) |
| Brio--Wu reference branch | не выбран | требуется решение владельца |
| публикация Бисикало--Жилкина | не идентифицирована | требуется точное название или PDF |
| cluster campaign | нет доступа/спецификации | требуется scheduler/account/partition/node data |
| шаблон НИРС | не предоставлен | draft использует нейтральный LaTeX layout |

## Допустимые утверждения после T00

- Реализованы CPU-путь AMReX, HLLD, MUSCL, staggered CT и SSP-RK2/Heun.
- Release-сборка и два standalone smoke-теста проходят на указанной машине.

## Запрещённые утверждения до следующих gate

- Глобальная консервативность AMR: нет gas reflux.
- Полная numerical verification: отсутствует CTest и независимые reference
  metrics.
- GPU port: горячие циклы используют `amrex::LoopOnCpu`.
- Scaling/ускорение: нет повторов, привязки ресурсов и raw timing protocol.
- Воспроизведение legacy/ВКР: input, mesh, output и exact run не связаны
  доказательным manifest.
