/**
 * @file rng_generators.hpp
 * @brief Три модифицированных генератора псевдослучайных чисел.
 * @details Содержит классы:
 *          - MLCG_XOR (модифицированный LCG с нелинейным перемешиванием)
 *          - LFG_XOR (генератор Фибоначчи с запаздыванием и XOR)
 *          - Xorshift128_MLT (модифицированный Xorshift128 с умножением)
 */

#pragma once

#include <cstdint>

/**
 * @brief Модифицированный LCG с нелинейным перемешиванием (MLCG-XOR).
 * @details Классический LCG: X_{n+1} = (a * X_n + c) mod 2^64.
 *          Модификации:
 *          - Модуль 2^64 за счёт естественного переполнения uint64_t.
 *          - Параметры a и c подобраны для полного периода 2^64.
 *          - После каждого шага применяется нелинейное перемешивание
 *            (XOR и сдвиги), устраняющее корреляцию младших битов.
 */
class MLCG_XOR {
public:
    using result_type = uint64_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789).
     */
    explicit MLCG_XOR(uint64_t seed = 123456789);

    /** @brief Генерирует 64-битное псевдослучайное число. */
    uint64_t next();

    /** @brief Минимальное значение, возвращаемое next(). */
    static constexpr uint64_t min() { return 0; }

    /** @brief Максимальное значение, возвращаемое next(). */
    static constexpr uint64_t max() { return ~0ULL; }

private:
    uint64_t state;                              ///< Текущее состояние
    static constexpr uint64_t A = 0x9e3779b97f4a7c15ULL; ///< Множитель (золотое сечение)
    static constexpr uint64_t C = 0xbf58476d1ce4e5b9ULL; ///< Приращение (нечётное)
};

/**
 * @brief Генератор Фибоначчи с запаздыванием и XOR (LFG-XOR).
 * @details Базовый алгоритм: LFG с лагами (55,24). Модификации:
 *          - Вместо сложения используется XOR (быстрее, лучше перемешивание).
 *          - Добавлена нелинейная функция scramble после обновлений.
 *          - Улучшенная инициализация буфера через LCG и перемешивание.
 */
class LFG_XOR {
public:
    using result_type = uint32_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789).
     */
    explicit LFG_XOR(uint32_t seed = 123456789);

    /** @brief Генерирует 32-битное псевдослучайное число. */
    uint32_t next();

    /** @brief Минимальное значение, возвращаемое next(). */
    static constexpr uint32_t min() { return 0; }

    /** @brief Максимальное значение, возвращаемое next(). */
    static constexpr uint32_t max() { return ~0U; }

private:
    static constexpr int LAG1 = 55;   ///< Первый лаг
    static constexpr int LAG2 = 24;   ///< Второй лаг
    uint32_t buffer[LAG1];            ///< Буфер состояний
    int index;                        ///< Текущий индекс

    /** @brief Нелинейное перемешивание 32-битного слова (XOR + сдвиги). */
    uint32_t scramble(uint32_t x) const;

    /** @brief Инициализация буфера с заданным seed. */
    void initBuffer(uint32_t seed);
};

/**
 * @brief Модифицированный Xorshift128 с нелинейным умножением (Xorshift128-MLT).
 * @details Базовый алгоритм: Xorshift128 (период 2^128-1). Модификации:
 *          - Добавлено нелинейное умножение на две константы
 *            и дополнительные XOR-сдвиги.
 *          - Прогрев (10 итераций) для устранения корреляции начальных значений.
 */
class Xorshift128_MLT {
public:
    using result_type = uint64_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789).
     */
    explicit Xorshift128_MLT(uint64_t seed = 123456789);

    /** @brief Генерирует 64-битное псевдослучайное число. */
    uint64_t next();

    /** @brief Минимальное значение, возвращаемое next(). */
    static constexpr uint64_t min() { return 0; }

    /** @brief Максимальное значение, возвращаемое next(). */
    static constexpr uint64_t max() { return ~0ULL; }

private:
    uint64_t s[2];   ///< Состояние (128 бит)
};
