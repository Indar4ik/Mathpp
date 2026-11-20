#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <print>
#include <string>
#include <utility>
#include <vector>

#include <raylib.h>

#include "../dft.hpp"

// Максимальная частота эпицикла
constexpr int M = 1024;
constexpr size_t NUM_EPICYCLES = 2 * M;
constexpr size_t AMORTIZATION_COUNT = 8;
constexpr double AMORTIZATION_DEEP_INV = 1.0 / static_cast<double>(1ULL << AMORTIZATION_COUNT);

// Тип комплексного числа
using Complex = std::complex<double>;

// Контур - последовательность точек
using Contour = std::vector<Complex>;

// Член ряда
struct Epicycle {
    double freq;  // Показатель экспоненты, характеризующий частоту вращения
    Complex coef; // Комплексный коэффициент (амплитуда и фаза)
    double radius;
};

// Отрисовка окружностей
struct Joint {
    Complex pos;
    double radius;
};

using EpicycleArray = std::array<Epicycle, NUM_EPICYCLES>;
using JointArray = std::array<Joint, NUM_EPICYCLES + 1>;

// Массив точек с изображения -> комплексный контур
Contour prepareContour(const std::vector<std::pair<int, int>>& raw_points) noexcept {
    if (raw_points.empty()) return {};

    Contour contour;
    contour.reserve(raw_points.size());

    for (const auto& pt : raw_points) contour.emplace_back(static_cast<double>(pt.first), static_cast<double>(pt.second));

    // Замыкание контура
    if (raw_points.back() != raw_points.front()) contour.emplace_back(static_cast<double>(raw_points.front().first), static_cast<double>(raw_points.front().second));

    return contour;
}

// Вставляет между каждой парой подряд идущих точек их 2^(AMORTIZATION_COUNT)-1 средних значений
void contourAmortization(Contour& contour) noexcept {
    if constexpr (AMORTIZATION_COUNT == 0) return;
    if (contour.empty()) return;

    const size_t n = contour.size();
    contour.resize(((contour.size() - 1) << AMORTIZATION_COUNT) + 1);

    contour.back() = contour[n - 1];
    for (size_t i = n - 1; i > 0; --i) {
        Complex x = contour[i - 1];
        Complex add = (contour[i] - x) * AMORTIZATION_DEEP_INV;
        size_t base = (i - 1) << AMORTIZATION_COUNT;

        #pragma clang loop unroll(full) vectorize(enable)
        for (size_t j = 0; j < (1 << AMORTIZATION_COUNT); ++j) {
            contour[base + j] = x;
            x += add;
        }
    }
}

// Центрирует замкнутый контур в (0, 0) и возвращает снятый центр.
// Последняя точка - дубликат первой (замыкание из prepareContour), в среднее она не входит:
// тогда среднее по Q различным вершинам амортизированной ломаной в точности равно
// c_0 = ∫z(t)dt - центру масс кривой, которую рисуют эпициклы.
Complex centerContour(Contour& contour) noexcept {
    const size_t N = contour.size();
    if (N == 0) return {};
    const size_t Q = N > 1 ? N - 1 : N;

    double sum_re = 0.0;
    double sum_im = 0.0;

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < Q; ++i) {
        sum_re += contour[i].real();
        sum_im += contour[i].imag();
    }

    double center_re = sum_re / static_cast<double>(Q);
    double center_im = sum_im / static_cast<double>(Q);

    // Сдвигаем все N точек, чтобы дубликат остался равен первой (resample_pow2 по нему узнаёт замыкание)
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < N; ++i) contour[i] = Complex(contour[i].real() - center_re, contour[i].imag() - center_im);

    return {center_re, center_im};
}

// Частоты ∊ [-M; M] \ 0
// Всего будет вычислено 2M эпициклов
EpicycleArray computeEpicycles(Contour& contour) noexcept {
    if (contour.size() < 2) return {}; // точка: все коэффициенты нулевые
    EpicycleArray epicycles;

    // N >= 2M+1, иначе диапазоны бинов [1; M] и [N-M; N) пересекаются: бин k попадёт и в частоту k, и в k-N.
    size_t N = 0;
    fft::resample_pow2(contour, N, 2 * M + 1);
    // std::println("{}", N);
    [[assume(N > 2 * M && (N & (N - 1ull)) == 0ull)]];
    const double N_inv = 1.0 / static_cast<double>(N);
    fft::fft(contour);

    for (size_t i = 1; i <= M; ++i) {
        epicycles[i - 1] = {static_cast<double>(i), contour[i] * N_inv, std::abs(contour[i]) * N_inv};
        // std::println("c_{}: {:.5f}", i, epicycles[i - 1].radius);
    }
    for (size_t i = N - M; i < N; ++i) {
        epicycles[2 * M + i - N] = {static_cast<double>(i) - static_cast<double>(N), contour[i] * N_inv, std::abs(contour[i]) * N_inv};
        // std::println("c_{}: {:.5f}", static_cast<int64_t>(i) - static_cast<int64_t>(N_new), epicycles[2 * M + i - N_new].radius);
    }

    std::ranges::sort(epicycles, std::ranges::greater{}, [](const Epicycle& c) { return c.radius; });

    return epicycles;
}

// Расчёт радиусов и позиций окружностей и векторов
void calculateFrame(EpicycleArray& epicycles, double t, JointArray& joints) noexcept {
    if (epicycles.empty()) return;

    const double pi2_t = 2.0 * std::numbers::pi * t;
    Complex cur_pos = {0.0, 0.0};

    // Без unroll(full): полная развёртка 2000 итераций давала 275 КБ кода (почти весь бинарник),
    // цикл не помещался ни в L1i, ни в µop-кэш и был в 1.32 раза медленнее.
    // Векторизовать его мешает префиксная сумма cur_pos (sin/cos clang векторизует только с -fveclib=libmvec).
    for (size_t i = 0; i < NUM_EPICYCLES; ++i) {
        double angle = pi2_t * epicycles[i].freq;

        double c = std::cos(angle);
        double s = std::sin(angle);

        double re = epicycles[i].coef.real();
        double im = epicycles[i].coef.imag();

        joints[i] = {cur_pos, epicycles[i].radius};
        cur_pos += Complex(re * c - im * s, re * s + im * c);
    }
    joints[NUM_EPICYCLES] = {cur_pos, 0.0};
}

std::vector<std::pair<int, int>> generateHeart() noexcept {
    std::vector<std::pair<int, int>> points;
    for (double t = 0; t < 2 * std::numbers::pi; t += 0.02) {
        // Формула сердечка
        double x = 240 * std::pow(std::sin(t), 3);
        double y = 15 * (13 * std::cos(t) - 5 * std::cos(2 * t) - 2 * std::cos(3 * t) - std::cos(4 * t));
        // std::print("({:.4f}, {:.4f}) ", x, y);
        points.emplace_back(static_cast<int>(x), static_cast<int>(y));
    }
    return points;
}
// const std::vector<std::pair<int, int>> raw_pixels = generateHeart();

std::vector<std::pair<int, int>> loadPoints(const std::string& filename) noexcept {
    std::vector<std::pair<int, int>> points;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error when open the file " << filename << "!\n";
        return points;
    }

    points.reserve(1000);

    int x, y;
    while (file >> x >> y) points.emplace_back(x, y);

    return points;
}

bool DrawToggleButton(Rectangle bounds, const char* text, bool& state) noexcept {
    Vector2 mouse = GetMousePosition();
    bool hovered = CheckCollisionPointRec(mouse, bounds);
    if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) state = !state;

    DrawRectangleRec(bounds, state ? (hovered ? DARKBLUE : BLUE) : (hovered ? LIGHTGRAY : GRAY));
    DrawRectangleLinesEx(bounds, 2, BLACK);
    int textWidth = MeasureText(text, 20);
    DrawText(text, bounds.x + (bounds.width - textWidth) / 2, bounds.y + 10, 20, WHITE);
    return hovered;
}

bool DrawClearButton(Rectangle bounds, const char* text, std::vector<Vector2>& vec) noexcept {
    Vector2 mouse = GetMousePosition();
    bool hovered = CheckCollisionPointRec(mouse, bounds);
    if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) vec.clear();

    DrawRectangleRec(bounds, hovered ? DARKBLUE : BLUE);
    DrawRectangleLinesEx(bounds, 2, BLACK);
    int textWidth = MeasureText(text, 20);
    DrawText(text, bounds.x + (bounds.width - textWidth) / 2, bounds.y + 10, 20, WHITE);
    return hovered;
}

// Реальная длина контура в пикселях: сумма длин всех отрезков исходного points.txt,
// включая замыкающий (последняя точка -> первая). Делает dt независимым от числа точек —
// скорость становится "пикселей контура за кадр", одинаковой для разреженных и плотных
// контуров (9 точек и 9000 точек идут с одинаковой видимой скоростью).
double contourPerimeter(const std::vector<std::pair<int, int>>& pts) noexcept {
    const size_t n = pts.size();
    if (n < 2) return 1.0;
    double len = 0.0;
    for (size_t i = 0; i < n - 1; ++i) {
        const auto& a = pts[i];
        const auto& b = pts[i + 1];
        const double dx = static_cast<double>(b.first - a.first);
        const double dy = static_cast<double>(b.second - a.second);
        len += std::sqrt(dx * dx + dy * dy);
    }
    const auto& a = pts.back();
    const auto& b = pts.front(); // замыкание контура
    const double dx = static_cast<double>(b.first - a.first);
    const double dy = static_cast<double>(b.second - a.second);
    len += std::sqrt(dx * dx + dy * dy);
    return len > 0.0 ? len : 1.0;
}

int main() {
    constexpr int width = 1920;
    constexpr int height = 1080;
    constexpr Vector2 offset = {width / 2.0f, height / 2.0f};

    const std::vector<std::pair<int, int>> raw_pixels = loadPoints("points.txt");
    if (raw_pixels.empty()) return 1;
    const double perimeter = contourPerimeter(raw_pixels);

    auto contour = prepareContour(raw_pixels);

    // Амортизация + центрирование
    contourAmortization(contour);
    const Complex center = centerContour(contour);

    // Силуэт сдвигается на тот же центр, что и кривая эпициклов, и замыкается как в prepareContour
    std::vector<Vector2> vec_pixels;
    vec_pixels.reserve(raw_pixels.size() + 1);
    for (const auto& pix : raw_pixels) vec_pixels.emplace_back(static_cast<float>(pix.first - center.real()) + offset.x, static_cast<float>(pix.second - center.imag()) + offset.y);
    vec_pixels.push_back(vec_pixels.front());

    auto epicycles = computeEpicycles(contour);

    // Преаллокация буферов для горячего цикла
    JointArray frame_joints;

    // Настраиваем антиалиасинг, чтобы линии были гладкими
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(width, height, "FourierPortrait");
    SetTargetFPS(240);

    double t = 0.0;
    double drawn = 0.0;                // накопленный параметр; 1.0 = один полный обход контура
    const double dt = 1.0 / perimeter; // шаг = 1 пиксель длины контура (speedMult -> px/кадр)
    std::vector<Vector2> trail;

    // --- ПЕРЕМЕННЫЕ ИНТЕРФЕЙСА ---
    bool showCircles = true;
    bool showLines = true;
    bool showTrail = true;
    bool showOutline = true;
    double speedMult = 1.0; // Начальная скорость

    while (!WindowShouldClose()) {
        // Управление скоростью стрелочками ВВЕРХ / ВНИЗ
        if (IsKeyPressed(KEY_UP)) speedMult += 0.05;
        if (IsKeyPressed(KEY_DOWN) && speedMult > 0.06) speedMult -= 0.05;
        if (IsKeyPressed(KEY_S)) speedMult = 0.004;

        // Физика
        calculateFrame(epicycles, t, frame_joints);

        double step = dt * speedMult;
        // Трасса пишется ровно один период независимо от скорости (даже если её менять по ходу).
        // Кнопка Clear опустошает trail -> накопитель сбрасывается и кривая перерисовывается заново.
        if (trail.empty()) drawn = 0.0;
        if (drawn < 1.0) {
            trail.emplace_back(static_cast<float>(frame_joints.back().pos.real()) + offset.x, static_cast<float>(frame_joints.back().pos.imag()) + offset.y);
            drawn += step;
        }

        t += step;
        t -= floor(t);

        // Отрисовка
        BeginDrawing();
        ClearBackground({20, 20, 20, 255}); // Темный фон

        // Рисуем эпициклы
        for (size_t i = 0; i < NUM_EPICYCLES; ++i) {
            Vector2 p1 = {static_cast<float>(frame_joints[i].pos.real()) + offset.x, static_cast<float>(frame_joints[i].pos.imag()) + offset.y};
            Vector2 p2 = {static_cast<float>(frame_joints[i + 1].pos.real()) + offset.x, static_cast<float>(frame_joints[i + 1].pos.imag()) + offset.y};
            float radius = static_cast<float>(frame_joints[i].radius);

            // Рисуем круги только если радиус > 0.5 пикселя, иначе их всё равно не видно
            if (showCircles && radius > 0.5f) DrawCircleLines(p1.x, p1.y, radius, Fade(WHITE, 0.12f));
            if (showLines) DrawLineV(p1, p2, WHITE);
        }

        // Рисуем след (polyline)
        if (showTrail && trail.size() > 1) DrawLineStrip(trail.data(), trail.size(), YELLOW);

        // Силуэт реальной кривой
        if (showOutline) DrawLineStrip(vec_pixels.data(), vec_pixels.size(), Fade(WHITE, 0.05f));

        // --- ОТРИСОВКА UI ---
        DrawToggleButton({20, 20, 140, 40}, "Circles", showCircles);
        DrawToggleButton({170, 20, 140, 40}, "Lines", showLines);
        DrawToggleButton({320, 20, 140, 40}, "Curve", showTrail);
        DrawToggleButton({470, 20, 140, 40}, "Outline", showOutline);
        DrawClearButton({620, 20, 140, 40}, "Clear", trail);

        // Текст со скоростью
        DrawText(TextFormat("Speed: %.2fx (Arrows UP/DOWN)", speedMult), 780, 30, 20, LIGHTGRAY);
        DrawText(TextFormat("Points: %zu | M: %d", trail.size(), M), 20, height - 30, 20, DARKGRAY);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
