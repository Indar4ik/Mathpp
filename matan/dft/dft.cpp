#include <SDL2/SDL_rect.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numbers>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>
#include "SDL2/SDL.h"

// Тип комплексного числа
using Complex = std::complex<double>;

// Контур - последовательность точек
using Contour = std::vector<Complex>;

// Член ряда
struct Epicycle {
    int freq; // Показатель экспоненты, характеризующий частоту вращения
    Complex coef; // Комплексный коэффициент (амплитуда и фаза)
    double radius;
};

// Отрисовка окружностей
struct Joint {
    Complex pos;
    double radius;
};

// Массив точек с изображения -> комплексный контур, центрированный на (0, 0)
Contour prepareContour(const std::vector<std::pair<double, double>>& raw_points) noexcept{
    if (raw_points.empty()) return {};

    const size_t n = raw_points.size();
    Contour contour;
    contour.reserve(n);
    for (const auto& pt : raw_points){
        contour.emplace_back(static_cast<double>(pt.first),
                             static_cast<double>(pt.second));
    }

    Complex sum = std::reduce(contour.begin(), contour.end(), Complex(0.0, 0.0));
    Complex center = sum / static_cast<double>(n);

    for (auto& z : contour){
        z -= center;
        // z = Complex(z.real(), -z.imag());
    }

    return contour;
}

// Частоты ∊ [-M; M]
// Всего будет вычислено 2M+1 эпициклов
std::vector<Epicycle> computeEpicycles(const Contour& contour, int M){
    const int N = contour.size();
    if (N == 0) return {};

    std::vector<Epicycle> epicycles;
    epicycles.reserve(2 * M + 1);

    const double pi2_N = -2.0 * std::numbers::pi / N;

    for (int n = -M; n <= M; ++n){
        Complex sum(0.0, 0.0);

        #pragma clang loop vectorize(enable)
        for (int k = 0; k < N; ++k){
            double angle = pi2_N * n * k;

            Complex w(std::cos(angle), std::sin(angle));

            sum += contour[k] * w;
        }

        epicycles.push_back({n, sum / static_cast<double>(N)});
    }

    std::ranges::sort(epicycles, std::ranges::greater{},[](const Epicycle& e) { return std::norm(e.coef); });

    return epicycles;
}

std::vector<Joint> calculateFrame(const std::vector<Epicycle>& epicycles, double t){
    if (epicycles.empty()) return {};

    const size_t n = epicycles.size();

    std::vector<Complex> vectors(n);

    const double pi2_t = 2.0 * std::numbers::pi * t;

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; ++i){
        const double angle = pi2_t * epicycles[i].freq;

        Complex rot(std::cos(angle), std::sin(angle));

        vectors[i] = rot * epicycles[i].coef;
    }

    std::vector<Joint> joints;
    joints.reserve(n + 1);

    Complex cur_pos(0.0, 0.0);

    for (size_t i = 0; i < n; ++i){
        joints.push_back({cur_pos, epicycles[i].radius});
        cur_pos += vectors[i];
    }

    joints.push_back({cur_pos, 0.0});

    return joints;
}

void DrawCircle(SDL_Renderer* renderer, float cx, float cy, float radius) {
    const int segments = 60;
    std::vector<SDL_FPoint> points(segments + 1);
    for (int i = 0; i <= segments; ++i) {
        float theta = 2.0f * std::numbers::pi * float(i) / float(segments);
        points[i] = { cx + radius * std::cos(theta), cy + radius * std::sin(theta) };
    }
    SDL_RenderDrawLinesF(renderer, points.data(), points.size());
}

std::vector<std::pair<double, double>> generateHeart() {
    std::vector<std::pair<double, double>> points;
    for (double t = 0; t < 2 * std::numbers::pi; t += 0.02) {
        // Формула сердечка
        double x = 16 * std::pow(std::sin(t), 3);
        double y = -(13 * std::cos(t) - 5 * std::cos(2*t) - 2 * std::cos(3*t) - std::cos(4*t));
        // Масштабируем
        points.push_back({x * 15, y * 10});
    }
    return points;
}

int main(){
    constexpr double speed = 0.1;
    // std::vector<std::pair<int, int>> raw_pixels = {{0, 0}, {10, 10}};
    auto raw_pixels = generateHeart();
    Contour contour = prepareContour(raw_pixels);

    constexpr int M = 100;
    auto epicycles = computeEpicycles(contour, M);

    const int width = 1280;
    const int height = 720;
    SDL_Window* window = SDL_CreateWindow("Fourier Epicycles (SDL2)", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    double t = 0.0;
    const double dt = 1.0 / contour.size();
    std::vector<SDL_FPoint> trail;
    float offsetX = width / 2.0f;
    float offsetY = height / 2.0f;
    bool running = true;

    // 3. Главный цикл
    while (running) {
        
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) running = false;
        }

        // Обновление физики
        auto joints = calculateFrame(epicycles, t);
        if (!joints.empty()) {
            trail.push_back({ (float)joints.back().pos.real() + offsetX, (float)joints.back().pos.imag() + offsetY });
        }

        t += dt * speed;
        t -= (int)t;

        // Отрисовка
        SDL_SetRenderDrawColor(renderer, 20, 20, 20, 255); // Темно-серый фон
        SDL_RenderClear(renderer);

        // Рисуем круги и линии (эпициклы)
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 200); // Полупрозрачный белый
        for (size_t i = 0; i < joints.size() - 1; ++i) {
            float cx = joints[i].pos.real() + offsetX;
            float cy = joints[i].pos.imag() + offsetY;
            float nx = joints[i+1].pos.real() + offsetX;
            float ny = joints[i+1].pos.imag() + offsetY;
            
            DrawCircle(renderer, cx, cy, joints[i].radius);
            SDL_RenderDrawLineF(renderer, cx, cy, nx, ny);
        }

        // Рисуем след карандаша
        if (trail.size() > 1) {
            SDL_SetRenderDrawColor(renderer, 255, 50, 50, 255); // Ярко-красный
            SDL_RenderDrawLinesF(renderer, trail.data(), trail.size());
        }

        SDL_RenderPresent(renderer); // Вывод кадра
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}