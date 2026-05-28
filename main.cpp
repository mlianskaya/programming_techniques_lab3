/**
 * @file main.cpp
 * @brief Программа для тестирования трёх модифицированных генераторов псевдослучайных чисел.
 
 * 1. Генерирует 20 выборок по 20 000 элементов для каждого генератора.
 * 2. Для каждой выборки вычисляет среднее, стандартное отклонение и коэффициент вариации.
 * 3. Проверяет равномерность распределения каждой выборки критерием хи-квадрат (уровень значимости 0.05).
 * 4. Применяет пять тестов NIST (частотный, блочный частотный, на серии, на длину серии, кумулятивных сумм).
 * 5. Замеряет время генерации от 1 000 до 1 000 000 чисел и сохраняет результаты в CSV.
 * 6. Результаты хи-квадрат и NIST выводятся в консоль.
 */

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <random>
#include "rng_generators.hpp"
#include "ChiSquared.hpp"
#include "NISTTests.hpp"

// Вспомогательные статистики

/**
 * @brief Вычисляет выборочное среднее.
 * @tparam T числовой тип элементов выборки.
 * @param sample вектор значений.
 * @return среднее арифметическое.
 */
template <typename T>
double mean(const std::vector<T>& sample) {
    return std::accumulate(sample.begin(), sample.end(), 0.0) / sample.size();
}

/**
 * @brief Вычисляет выборочное стандартное отклонение.
 * @tparam T числовой тип элементов выборки.
 * @param sample вектор значений.
 * @return стандартное отклонение (корень из дисперсии).
 */
template <typename T>
double stddev(const std::vector<T>& sample) {
    double m = mean(sample);
    double sq = 0.0;
    for (const auto& x : sample) sq += (x - m) * (x - m);
    return std::sqrt(sq / sample.size());
}

/**
 * @brief Вычисляет коэффициент вариации (относительное стандартное отклонение).
 * @tparam T числовой тип элементов выборки.
 * @param sample вектор значений.
 * @return отношение стандартного отклонения к среднему.
 */
template <typename T>
double variability(const std::vector<T>& sample) {
    return stddev(sample) / mean(sample);
}

// Преобразование выборки в битовую строку 

/**
 * @brief Преобразует вектор целых чисел в вектор битов (младшие биты сначала).
 * @tparam T целочисленный тип (должен поддерживать побитовые операции).
 * @param sample вектор чисел.
 * @param bitsPerNumber сколько младших битов использовать (обычно sizeof(T)*8).
 * @return вектор bool, представляющий битовую строку.
 */
template <typename T>
std::vector<bool> toBits(const std::vector<T>& sample, int bitsPerNumber) {
    std::vector<bool> bits;
    bits.reserve(sample.size() * bitsPerNumber);
    for (T x : sample) {
        for (int b = 0; b < bitsPerNumber; ++b) {
            bits.push_back((x >> b) & 1);
        }
    }
    return bits;
}

// Генерация выборок

/**
 * @brief Генерирует заданное количество выборок указанного размера с помощью переданного генератора.
 * @tparam Generator тип генератора (должен иметь метод next()).
 * @param gen ссылка на генератор.
 * @param numSamples количество выборок.
 * @param sampleSize размер каждой выборки.
 * @return вектор векторов сгенерированных чисел.
 */
template <typename Generator>
std::vector<std::vector<typename Generator::result_type>>
generateSamples(Generator& gen, size_t numSamples, size_t sampleSize) {
    std::vector<std::vector<typename Generator::result_type>> samples(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i].reserve(sampleSize);
        for (size_t j = 0; j < sampleSize; ++j) {
            samples[i].push_back(gen.next());
        }
    }
    return samples;
}

// Вывод статистики и хи-квадрат 

/**
 * @brief Выводит для каждой выборки среднее, стандартное отклонение, коэффициент вариации
 *        и результат критерия хи-квадрат.
 * @tparam Generator тип генератора (используется для получения min/max).
 * @param name имя генератора.
 * @param samples вектор выборок.
 * @param significance уровень значимости для хи-квадрат теста.
 */
template <typename T>
void printStats(const std::string& name,
                const std::vector<std::vector<T>>& samples,
                T minVal, T maxVal,
                double significance) {
    std::cout << "\n=== " << name << " ===\n";
    size_t passed = 0;
    for (size_t i = 0; i < samples.size(); ++i) {
        double m = mean(samples[i]);
        double sd = stddev(samples[i]);
        double cv = variability(samples[i]);
        bool chi_ok = custom::chi_squared_test(samples[i],
                                              minVal, maxVal,
                                              significance,
                                              custom::DistributionType::UNIFORM_INT);
        if (chi_ok) ++passed;
        std::cout << "Sample " << std::setw(2) << i+1
                  << ": mean=" << std::fixed << std::setprecision(2) << m
                  << ", std=" << sd
                  << ", var.coef=" << cv
                  << ", Chi-square=" << (chi_ok ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "Chi-square passed: " << passed << "/" << samples.size() << "\n";
}

// NIST тесты 

/**
 * @brief Применяет пять тестов NIST ко всем выборкам (объединяя их в одну битовую последовательность).
 * @tparam Generator тип генератора (используется для получения размера типа).
 * @param name имя генератора.
 * @param samples вектор выборок.
 * @param significance уровень значимости для NIST тестов.
 */
template <typename T>
void runNISTTests(const std::string& name,
                  const std::vector<std::vector<T>>& samples,
                  double significance) {
    std::cout << "\n--- NIST tests for " << name << " ---\n";
    std::vector<bool> allBits;
    for (const auto& sample : samples) {
        auto bits = toBits(sample, sizeof(T) * 8);
        allBits.insert(allBits.end(), bits.begin(), bits.end());
    }

    nist::FrequencyTest freq(allBits);
    nist::BlockFrequencyTest blockFreq(allBits, 8);
    nist::RunsTest runs(allBits);
    nist::LongestRunTest longestRun(allBits);
    nist::CumulativeSumsTest cumSums(allBits);

    auto status = [&](double p) { return (p >= significance) ? "PASS" : "FAIL"; };
    std::cout << "Frequency test:          " << status(freq.p_value) << "\n";
    std::cout << "Block frequency test:    " << status(blockFreq.p_value) << "\n";
    std::cout << "Runs test:               " << status(runs.p_value) << "\n";
    std::cout << "Longest run test:        " << status(longestRun.p_value) << "\n";
    std::cout << "Cumulative sums test:    " << status(cumSums.p_value) << "\n";
}

// Замеры времени 

/**
 * @brief Замеряет время генерации заданного количества чисел для переданного генератора.
 * @tparam Generator тип генератора.
 * @param gen ссылка на генератор.
 * @param sizes вектор размеров (количеств генерируемых чисел).
 * @return вектор времени в наносекундах для каждого размера.
 */
template <typename Generator>
std::vector<uint64_t> measureTimes(Generator& gen, const std::vector<size_t>& sizes) {
    std::vector<uint64_t> times;
    for (size_t n : sizes) {
        auto start = std::chrono::high_resolution_clock::now();
        volatile uint64_t tmp;       // volatile предотвращает оптимизацию
        for (size_t j = 0; j < n; ++j) tmp = gen.next();
        auto end = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    }
    return times;
}


int main() {
    // Параметры эксперимента (соответствуют заданию)
    const size_t NUM_SAMPLES = 20;        // не менее 20 выборок
    const size_t SAMPLE_SIZE = 20000;     // каждая выборка ≥1000 элементов
    const double SIGNIFICANCE = 0.05;     // уровень значимости 0.05

    // 1. Три модифицированных генератора 
    MLCG_XOR       g1(123456789);
    LFG_XOR        g2(123456789);
    Xorshift128_MLT g3(123456789);

    // 2. Генерация выборок 
    auto samples1 = generateSamples(g1, NUM_SAMPLES, SAMPLE_SIZE);
    auto samples2 = generateSamples(g2, NUM_SAMPLES, SAMPLE_SIZE);
    auto samples3 = generateSamples(g3, NUM_SAMPLES, SAMPLE_SIZE);

    // 3. Статистика и хи-квадрат для каждой выборки 
    printStats("MLCG-XOR", samples1, MLCG_XOR::min(), MLCG_XOR::max(), SIGNIFICANCE);
    printStats("LFG-XOR", samples2, LFG_XOR::min(), LFG_XOR::max(), SIGNIFICANCE);
    printStats("Xorshift128-MLT", samples3, Xorshift128_MLT::min(), Xorshift128_MLT::max(), SIGNIFICANCE);

    // 4. NIST тесты (пять тестов) 
    runNISTTests("MLCG-XOR", samples1, SIGNIFICANCE);
    runNISTTests("LFG-XOR", samples2, SIGNIFICANCE);
    runNISTTests("Xorshift128-MLT", samples3, SIGNIFICANCE);

    // 5. Замеры времени генерации 
    std::vector<size_t> sizes = {1000, 50000, 100000, 500000, 1000000};
    // Создаём новые экземпляры генераторов
    MLCG_XOR       t1(123456789);
    LFG_XOR        t2(123456789);
    Xorshift128_MLT t3(123456789);

    auto timings1 = measureTimes(t1, sizes);
    auto timings2 = measureTimes(t2, sizes);
    auto timings3 = measureTimes(t3, sizes);

    // Сохраняем результаты в CSV для построения графиков
    std::ofstream csv("timings.csv");
    csv << "size,MLCG_XOR,LFG_XOR,Xorshift128_MLT\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        csv << sizes[i] << ","
            << timings1[i] << ","
            << timings2[i] << ","
            << timings3[i] << "\n";
    }
    csv.close();

    std::cout << "\nTiming results saved to timings.csv\n";
    return 0;
}