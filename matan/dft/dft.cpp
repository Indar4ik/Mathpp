#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <numeric>
// #include <ranges>
// #include <string>
#include <utility>
#include <vector>
#include "raylib/include/raylib.h"

// Тип комплексного числа
using Complex = std::complex<double>;

// Контур - последовательность точек
using Contour = std::vector<Complex>;

// Член ряда
struct Epicycle {
    double freq; // Показатель экспоненты, характеризующий частоту вращения
    Complex coef; // Комплексный коэффициент (амплитуда и фаза)
    double radius;
};

// Отрисовка окружностей
struct Joint {
    Complex pos;
    double radius;
};

// Массив точек с изображения -> комплексный контур, центрированный на (0, 0)
Contour prepareContour(const std::vector<std::pair<int, int>>& raw_points) noexcept{
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
        z = Complex(z.real(), -z.imag());
    }

    return contour;
}

// Частоты ∊ [-M; M]
// Всего будет вычислено 2M+1 эпициклов
std::vector<Epicycle> computeEpicycles(const Contour& contour, int M) noexcept{
    const int N = contour.size();
    if (N == 0) return {};

    std::vector<Epicycle> epicycles;
    epicycles.reserve(2 * M + 1);

    const double pi2_N = -2.0 * std::numbers::pi / N;

    for (int n = -M; n <= M; ++n){
        Complex sum(0.0, 0.0);

        #pragma clang loop vectorize(enable)
        for (size_t k = 0; k < N; ++k){
            double angle = pi2_N * n * k;
            sum += contour[k] * Complex(std::cos(angle), std::sin(angle));
        }
        
        Complex c_n = sum / static_cast<double>(N);
        epicycles.push_back({static_cast<double>(n), c_n, std::abs(c_n)});
    }

    std::ranges::sort(epicycles, std::ranges::greater{},[](const Epicycle& e) { return std::norm(e.coef); });

    return epicycles;
}

std::vector<Joint> calculateFrame(const std::vector<Epicycle>& epicycles, double t) noexcept{
    if (epicycles.empty()) return {};

    const size_t n = epicycles.size();

    std::vector<Complex> vectors(n);

    const double pi2_t = 2.0 * std::numbers::pi * t;

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; ++i){
        double angle = pi2_t * epicycles[i].freq;
        
        double c = std::cos(angle);
        double s = std::sin(angle);

        double re = epicycles[i].coef.real();
        double im = epicycles[i].coef.imag();

        vectors[i] = Complex(re * c - im * s, re * s + im * c);
    }

    std::vector<Joint> joints(n + 1);

    Complex cur_pos = {0.0, 0.0};
    for (size_t i = 0; i < n; ++i){
        joints[i] = {cur_pos, epicycles[i].radius};
        cur_pos += vectors[i];
    }
    joints[n] = {cur_pos, 0.0};

    return joints;
}

std::vector<std::pair<int, int>> generateHeart() {
    std::vector<std::pair<int, int>> points;
    for (double t = 0; t < 2 * std::numbers::pi; t += 0.02) {
        // Формула сердечка
        double x = 240 * std::pow(std::sin(t), 3);
        double y = 15 * (13 * std::cos(t) - 5 * std::cos(2*t) - 2 * std::cos(3*t) - std::cos(4*t));
        // Масштабируем
        points.push_back({(int)x, (int)y});
    }
    return points;
}

bool DrawToggleButton(Rectangle bounds, const char* text, bool& state) {
    Vector2 mouse = GetMousePosition();
    bool hovered = CheckCollisionPointRec(mouse, bounds);
    if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) state = !state;
    
    DrawRectangleRec(bounds, state ? (hovered ? DARKBLUE : BLUE) : (hovered ? LIGHTGRAY : GRAY));
    DrawRectangleLinesEx(bounds, 2, BLACK);
    int textWidth = MeasureText(text, 20);
    DrawText(text, bounds.x + (bounds.width - textWidth)/2, bounds.y + 10, 20, WHITE);
    return hovered;
}

int main(){
    constexpr int M = 5;
    constexpr int width = 1280;
    constexpr int height = 720;
    constexpr Vector2 offset = { width / 2.0f, height / 2.0f };

    const std::vector<std::pair<int, int>> raw_pixels = generateHeart();
    const double pixels_count = raw_pixels.size();
    auto contour = prepareContour(raw_pixels);
    std::vector<Vector2> vec_pixels;
    vec_pixels.reserve(raw_pixels.size());
    for (auto& pix : contour){
        vec_pixels.push_back(Vector2((float)pix.real() + offset.x, (float)pix.imag() + offset.y));
    }
    auto epicycles = computeEpicycles(contour, M);
    
    // Настраиваем антиалиасинг, чтобы линии были гладкими
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(width, height, "Raylib Fourier Epicycles");
    SetTargetFPS(240);

    double t = 0.0;
    const double dt = 1.0 / pixels_count;
    std::vector<Vector2> trail;

    // --- ПЕРЕМЕННЫЕ ИНТЕРФЕЙСА ---
    bool showCircles = true;
    bool showLines = true;
    bool showTrail = true;
    double speedMult = 0.25; // Начальная скорость

    while (!WindowShouldClose()) {
        // Управление скоростью стрелочками ВВЕРХ / ВНИЗ
        if (IsKeyPressed(KEY_UP)) speedMult += 0.05;
        if (IsKeyPressed(KEY_DOWN) && speedMult > 0.05) speedMult -= 0.05;

        // Физика
        auto joints = calculateFrame(epicycles, t);
        if (!joints.empty() && showTrail && (speedMult * trail.size() < pixels_count)) {
            trail.push_back({ (float)joints.back().pos.real() + offset.x, (float)joints.back().pos.imag() + offset.y });
        }

        t += dt * speedMult;
        t -= (int)t;

        // Отрисовка
        BeginDrawing();
        ClearBackground({20, 20, 20, 255}); // Темный фон

        // Рисуем эпициклы
        for (size_t i = 0; i < joints.size() - 1; ++i) {
            Vector2 p1 = { (float)joints[i].pos.real() + offset.x, (float)joints[i].pos.imag() + offset.y };
            Vector2 p2 = { (float)joints[i+1].pos.real() + offset.x, (float)joints[i+1].pos.imag() + offset.y };
            float radius = (float)joints[i].radius;

            // Рисуем круги только если радиус > 0.5 пикселя, иначе их всё равно не видно
            if (showCircles && radius > 0.5f) {
                DrawCircleLines(p1.x, p1.y, radius, Fade(WHITE, 0.2f));
            }
            if (showLines) {
                DrawLineV(p1, p2, Fade(WHITE, 0.6f));
            }
        }

        // Рисуем след (polyline)
        if (showTrail && trail.size() > 1) {
            DrawLineStrip(trail.data(), trail.size(), RED);
        }

        // Силуэт реальной кривой
        DrawLineStrip(vec_pixels.data(), vec_pixels.size(), {0, 0, 255, 100});

        // --- ОТРИСОВКА UI ---
        DrawToggleButton({20, 20, 140, 40}, "Circles", showCircles);
        DrawToggleButton({170, 20, 140, 40}, "Lines", showLines);
        DrawToggleButton({320, 20, 140, 40}, "Curve", showTrail);
        
        // Текст со скоростью
        DrawText(TextFormat("Speed: %.2fx (Arrows UP/DOWN)", speedMult), 480, 30, 20, LIGHTGRAY);
        DrawText(TextFormat("Points: %zu | M: %d", trail.size(), M), 20, height - 30, 20, DARKGRAY);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}