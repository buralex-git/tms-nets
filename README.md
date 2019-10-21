СПЕЦИФИКАЦИЯ ПРОЕКТА
====================

## Описание

Данный репозиторий содержит в себе программу генерации (t, m, s)-сетей с основанием 2, реализованную с помощью алгоритма,
предложенного Гаральдом Нидеррайтером в 1987 году в работе "Low-Discrepancy and Low-Dispersion Sequences".

## Ветви репозитория

  * **master** - главная ветвь; в ней располагается последняя стабильная версия проекта;
  * **development** - ветвь, содержащая самые актуальные наработки;
  * **gh-pages** - техническая ветвь, хранящая в себе файлы документации;
  * **nets-tests** - ветвь, содержащая модифицированный код генерации (t, m, s)-сетей и функции тестирования.

Объединение **gh-pages** с другими ветвями, а также самостоятельное изменение её содержимого не допускается.
Объединение **nets-tests** с другими ветвями и использование исходного кода (в особенности, методов генерации сетей) из этой
ветви в других ветвях допускается _только_ при полном понимании последствий производимых действий.

## Документация к исходному коду

Документация к стабильной версии программы доступна по адресу: https://jointpoints.github.io/tms-nets/ (на английском языке),
генерируется автоматически с помощью утилиты Doxygen при каждом изменении в ветви **master**.

## Компиляция и сборка

Для компиляции программы потребуются папки include и src. Репозиторий включает в себя только исходные файлы, которые могут быть
импортированы в любую IDE.

Программа стабильно компилируется на gcc; известны незначительные трудности с Visual C++, связанные с функцией `timestamp`.