/**
 * @file ChiSquared.hpp
 * @brief Реализация критерия согласия хи-квадрат (Пирсона) для проверки равномерности распределения.
 * @details Содержит функции chi_squared (вычисление статистики и степеней свободы)
 *          и chi_squared_test (проверка гипотезы на заданном уровне значимости).
 */

#pragma once

#include <vector>
#include <algorithm>   
#include <cmath>       
#include <stdexcept>   
#include <sstream>     

namespace custom {

/**
 * @brief Тип распределения: целочисленное равномерное или непрерывное равномерное.
 */

enum class DistributionType { 
    UNIFORM_INT,   ///< Дискретное равномерное на [min, max] (целые)
    UNIFORM_REAL   ///< Непрерывное равномерное на [min, max) (вещественные)
};

/**
 * @brief Вычисляет статистику хи-квадрат и число степеней свободы для выборки.
 * @tparam Number Числовой тип элементов выборки (должен поддерживать сравнение, преобразование к double).
 * @param sample Вектор наблюдаемых значений.
 * @param min Минимальное возможное значение (включительно).
 * @param max Максимальное возможное значение (включительно для UNIFORM_INT, исключительно для UNIFORM_REAL).
 * @param type Тип распределения (целочисленное или непрерывное).
 * @return std::pair<double, std::size_t> Первый элемент — значение хи-квадрат, второй — число степеней свободы.
 * @throws std::invalid_argument если min >= max, или размер выборки < 2, или значение выборки вне диапазона.
 */

template <typename Number>
std::pair<double, std::size_t> chi_squared(const std::vector<Number>& sample,
                                           Number min, Number max,
                                           DistributionType type) {
    // Проверка корректности входных параметров
    if (min >= max) {
        std::ostringstream oss;
        oss << "Range minimum " << min << " is not less than range maximum " << max;
        throw std::invalid_argument(oss.str());
    }
    if (sample.size() < 2) {
        std::ostringstream oss;
        oss << "Sample size " << sample.size() << " is too small (need at least 2)";
        throw std::invalid_argument(oss.str());
    }

    //Определяем длину интервала
    int add = (type == DistributionType::UNIFORM_INT) ? 1 : 0;
    double interval = static_cast<double>(max) - min + add;

    // Количество интервалов (карманов) k = floor(log2(n)) + 1
    std::size_t num_buckets = static_cast<size_t>(std::log2(sample.size()) + 1);
    // Для целочисленных: если интервал меньше числа карманов, ограничиваем числом возможных значений
    if (type == DistributionType::UNIFORM_INT && static_cast<double>(num_buckets) > interval)
        num_buckets = static_cast<size_t>(interval);
    
    double bucket_size = interval / num_buckets;   // ширина одного интервала

    // Вычисляем границы интервалов
    std::vector<double> edges(num_buckets + 1);
    edges[0] = min;
    for (size_t i = 1; i <= num_buckets; ++i)
        edges[i] = min + bucket_size * i;

    // Теоретические вероятности для каждого интервала (равномерное распределение) 
    // = (ширина интервала)/interval
    std::vector<double> p_unif(num_buckets);
    for (size_t i = 0; i < num_buckets; ++i)
        p_unif[i] = (edges[i+1] - edges[i]) / interval;

    // Наблюдаемые частоты: подсчёт количества элементов в каждом интервале
    // obs[i] — количество элементов выборки, попавших в i-й интервал
    std::vector<size_t> obs(num_buckets, 0);
    for (Number val : sample) {
        // Ищем позицию val в массиве границ (первый элемент, строго больший val)
        auto it = std::upper_bound(edges.begin(), edges.end(), val);
        size_t pos = std::distance(edges.begin(), it) - 1;
        // Защита от выхода за границы (если val точно равно правой границе последнего интервала)
        if (pos >= num_buckets) pos = num_buckets - 1;
        ++obs[pos];
    }

    // Вычисление статистики хи-квадрат
    // χ² = Σ (obs_i - exp_i)^2 / exp_i
    // Возвращаем χ² и степени свободы df = k - 1
    double chi = 0.0;
    for (size_t i = 0; i < num_buckets; ++i) {
        double expected = p_unif[i] * sample.size();
        double diff = obs[i] - expected;
        chi += diff * diff / expected;
    }
    return {chi, num_buckets - 1};
}

/**
 * @brief Проверяет гипотезу о равномерном распределении выборки по критерию хи-квадрат.
 * @tparam Number Тип элементов выборки.
 * @param sample Вектор наблюдаемых значений.
 * @param min Минимальное возможное значение.
 * @param max Максимальное возможное значение.
 * @param significance Уровень значимости. Должен быть в интервале (0,1).
 * @param type Тип распределения (UNIFORM_INT или UNIFORM_REAL).
 * @return true, если гипотеза о равномерном распределении не отвергается на данном уровне значимости.
 */
template <typename Number>
bool chi_squared_test(const std::vector<Number>& sample,
                      Number min, Number max, double significance,
                      DistributionType type) {
    // Уровень значимости должен быть строго между 0 и 1
    if (significance <= 0 || significance >= 1) {
        std::ostringstream oss;
        oss << "Significance level " << significance << " is not in (0; 1) interval";
        throw std::invalid_argument(oss.str());
    }

    //Вычисляем χ² и df с помощью chi_squared
    auto [chi, df] = chi_squared(sample, min, max, type);

    // Для больших df используем нормальную аппроксимацию (центральная предельная теорема)
    // χ² приблизительно нормально с параметрами df, 2df
    // Возвращаем true, если p-value > significance (или если χ² меньше критического значения)
    if (df > 30) {
        // Нормированное отклонение: (χ² - df) / √(2df)
        double z = (chi - df) / std::sqrt(2 * df);
        // p-value = 1 - Φ(z) = 0.5 * erfc(z / √2)
        // Ф — функция распределения стандартного нормального закона
        double p = 0.5 * std::erfc(z / std::sqrt(2.0));
        return p > significance;
    } else {
        // Для малых df используем заранее вычисленные критические значения для уровня 0.05
        static const double critical[] = {3.84, 5.99, 7.81, 9.49, 11.07, 12.59, 14.07, 15.51, 16.92, 18.31};
        if (df <= 10)
            return chi <= critical[df - 1];
        else
            // Приближённая формула для df > 10 (до 30): χ² ≤ df + 2√(2df)
            return chi <= df + 2 * std::sqrt(2 * df);
    }
}

} 
