# T00 decision register

Дата: 2026-08-26. `owner` означает сторону, которая должна дать внешний
вход; не является назначением человека без его согласия.

| ID | Решение или блокер | Текущее состояние | Owner / required input | Разблокирует |
|---|---|---|---|---|
| D-001 | Официальный historical legacy commit | кандидат `9d0f60e…`, не подтверждён владельцем | owner: подтверждение commit/tag | T03 |
| D-002 | Provenance legacy Brio--Wu result | отсутствует config/run manifest | owner: exact input/config либо разрешение на новый run | T03, T05, T08 |
| D-003 | Brio--Wu reference branch | не выбрана | owner: regular или compound/non-regular reference и источник | T06 |
| D-004 | Публикация Бисикало--Жилкина | не идентифицирована | owner: название, DOI или PDF | report literature |
| D-005 | Кластер | нет доступа и параметров | owner: scheduler, account, partition, limits, node/GPU specs | T13 |
| D-006 | Формат НИРС | шаблон не предоставлен | owner: template/format rules, если обязательны | final T14 formatting |
| D-007 | AMR conservation | известный дефект: нет gas flux reflux | internal: реализовать и проверить отдельной фазой | T07 |
| D-008 | GPU status | host-only loops | internal: complete CUDA work package and parity | T12 |

## Зафиксированные решения

1. Первые три фазы означают T00, T01 и T02; T03 не начинается без D-001/D-002.
2. Текущая временная формула остаётся SSP-RK2/Heun (TVD RK2 Шу--Ошера), а не
   midpoint RK2.
3. T01 описывает текущие результаты как `implemented` или `diagnostic`; не
   повышает их статус до `verified` или `measured benchmark`.
4. T02 не меняет HLLD/MUSCL/CT и не создаёт performance claims.
