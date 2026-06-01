#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
// #include <mkl.h>
#include <numbers>
// #include <print>
#include <string>
#include <utility>
#include <vector>
#include <raylib.h>

// Максимальная частота эпицикла
constexpr int M = 400;
constexpr size_t NUM_EPICYCLES = 2 * M + 1;
constexpr size_t AMORTIZATION_COUNT = 10;
constexpr double AMORTIZATION_DEEP_INV = 1.0 / static_cast<double>(1ULL << AMORTIZATION_COUNT);

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

using EpicycleArray = std::array<Epicycle, NUM_EPICYCLES>;
using VectorArray = std::array<Complex, NUM_EPICYCLES>;
using JointArray = std::array<Joint, NUM_EPICYCLES + 1>;

// Массив точек с изображения -> комплексный контур
Contour prepareContour(const std::vector<std::pair<int, int>>& raw_points) noexcept {
    if (raw_points.empty()) return {};

    Contour contour;
    contour.reserve(raw_points.size());

    for (const auto& pt : raw_points){
        contour.emplace_back(static_cast<double>(pt.first), static_cast<double>(pt.second));
    }

    // Замыкание контура
    if(raw_points.back() != raw_points.front()) contour.emplace_back(static_cast<double>(raw_points.front().first), static_cast<double>(raw_points.front().second));

    return contour;
}

// Вставляет между каждой парой подряд идущих точек их 2^(AMORTIZATION_COUNT)-1 средних значений
void contourAmortization(Contour& contour) noexcept {
    if (contour.empty()) return;

    const size_t n = contour.size();
    contour.resize(((contour.size() - 1) << AMORTIZATION_COUNT) + 1);
    
    contour.back() = contour[n - 1];
    for (size_t i = n - 1; i > 0; --i){
        Complex x = contour[i - 1];
        const Complex add = (contour[i] - x) * AMORTIZATION_DEEP_INV;

        #pragma clang loop unroll(full) vectorize(enable)
        for (size_t j = 0; j < (1 << AMORTIZATION_COUNT); ++j){
            contour[((i - 1) << AMORTIZATION_COUNT) + j] = x;
            x += add;
        }
    }
}

// Центрирует контур в (0, 0)
void centerContour(Contour& contour) noexcept {
    const size_t N = contour.size();
    if (N == 0) return;

    double sum_re = 0.0;
    double sum_im = 0.0;

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < N; ++i) {
        sum_re += contour[i].real();
        sum_im += contour[i].imag();
    }

    double center_re = sum_re / static_cast<double>(N);
    double center_im = sum_im / static_cast<double>(N);

    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < N; ++i) {
        contour[i] = Complex(contour[i].real() - center_re, contour[i].imag() - center_im);
    }
}

// Частоты ∊ [-M; M]
// Всего будет вычислено 2M+1 эпициклов
EpicycleArray computeEpicycles(const Contour& contour) noexcept {
    const size_t N = contour.size();
    EpicycleArray epicycles;
    if (N == 0) return {};
    epicycles[M] = {0, 0, 0};

    const double pi2_N = -2.0 * std::numbers::pi / N;

    #pragma clang loop unroll(full)
    for (int n = 1; n <= M; ++n) {
        double sum_re1 = 0.0;
        double sum_im1 = 0.0;
        double sum_re2 = 0.0;
        double sum_im2 = 0.0;

        #pragma clang loop vectorize(enable)
        for (size_t k = 0; k < N; ++k) {
            double angle = pi2_N * n * static_cast<double>(k);
            
            double c = std::cos(angle);
            double s = std::sin(angle);
            
            double re = contour[k].real();
            double im = contour[k].imag();

            sum_re1 += re * c - im * s;
            sum_im1 += re * s + im * c;
            sum_re2 += re * c + im * s;
            sum_im2 += im * c - re * s;
        }
        
        Complex c_n1(sum_re1 / N, sum_im1 / N);
        Complex c_n2(sum_re2 / N, sum_im2 / N);
        // std::println("c_{}: {:.5f}", n, std::abs(c_n1));
        // std::println("c_{}: {:.5f}", -n, std::abs(c_n2));
        epicycles[M + n] = {static_cast<double>(n), c_n1, std::abs(c_n1)};
        epicycles[M - n] = {static_cast<double>(-n), c_n2, std::abs(c_n2)};
    }

    std::ranges::sort(epicycles, std::ranges::greater{},[](const Epicycle& e) { return std::norm(e.coef); });

    return epicycles;
}

// Расчёт радиусов и позиций окружностей и векторов
void calculateFrame(EpicycleArray& epicycles, double t, VectorArray& vectors, JointArray& joints) noexcept {
    if (epicycles.empty()) return;

    const double pi2_t = 2.0 * std::numbers::pi * t;

    #pragma clang loop vectorize(enable) unroll(full)
    for (size_t i = 0; i < NUM_EPICYCLES; ++i){
        double angle = pi2_t * epicycles[i].freq;
        
        double c = std::cos(angle);
        double s = std::sin(angle);

        double re = epicycles[i].coef.real();
        double im = epicycles[i].coef.imag();

        vectors[i] = Complex(re * c - im * s, re * s + im * c);
    }

    Complex cur_pos = {0.0, 0.0};

    #pragma clang loop unroll(full)
    for (size_t i = 0; i < NUM_EPICYCLES; ++i){
        joints[i] = {cur_pos, epicycles[i].radius};
        cur_pos += vectors[i];
    }
    joints[NUM_EPICYCLES] = {cur_pos, 0.0};
}

std::vector<std::pair<int, int>> generateHeart() noexcept {
    std::vector<std::pair<int, int>> points;
    for (double t = 0; t < 2 * std::numbers::pi; t += 0.02) {
        // Формула сердечка
        double x = 240 * std::pow(std::sin(t), 3);
        double y = 15 * (13 * std::cos(t) - 5 * std::cos(2*t) - 2 * std::cos(3*t) - std::cos(4*t));
        // std::print("({:.4f}, {:.4f}) ", x, y);
        points.push_back({(int)x, (int)y});
    }
    return points;
}

std::vector<std::pair<int, int>> loadPoints(const std::string& filename) noexcept {
    std::vector<std::pair<int, int>> points;
    std::ifstream file(filename);

    if (!file.is_open()){
        std::cerr << "Error when open the file " << filename << "!\n";
        return points;
    }

    points.reserve(1000);

    int x, y;
    while (file >> x >> y){
        // std::print("({}, {}) ", x, y);
        points.emplace_back(x, y);
    }

    return points;
}

bool DrawToggleButton(Rectangle bounds, const char* text, bool& state) noexcept {
    Vector2 mouse = GetMousePosition();
    bool hovered = CheckCollisionPointRec(mouse, bounds);
    if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) state = !state;
    
    DrawRectangleRec(bounds, state ? (hovered ? DARKBLUE : BLUE) : (hovered ? LIGHTGRAY : GRAY));
    DrawRectangleLinesEx(bounds, 2, BLACK);
    int textWidth = MeasureText(text, 20);
    DrawText(text, bounds.x + (bounds.width - textWidth)/2, bounds.y + 10, 20, WHITE);
    return hovered;
}

bool DrawClearButton(Rectangle bounds, const char* text, std::vector<Vector2>& vec) noexcept {
    Vector2 mouse = GetMousePosition();
    bool hovered = CheckCollisionPointRec(mouse, bounds);
    if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) vec.clear();
    
    DrawRectangleRec(bounds, hovered ? DARKBLUE : BLUE);
    DrawRectangleLinesEx(bounds, 2, BLACK);
    int textWidth = MeasureText(text, 20);
    DrawText(text, bounds.x + (bounds.width - textWidth)/2, bounds.y + 10, 20, WHITE);
    return hovered;
}

int main(){
    constexpr int width = 1280;
    constexpr int height = 720;
    constexpr Vector2 offset = { width / 2.0f, height / 2.0f };

    const std::vector<std::pair<int, int>> raw_pixels = loadPoints("points.txt");
    // const std::vector<std::pair<int, int>> raw_pixels = generateHeart();
    if (raw_pixels.empty()) return 1;
    const double pixels_count = static_cast<double>(raw_pixels.size());

    auto contour = prepareContour(raw_pixels);

    // Амортизация + центрирование
    if (AMORTIZATION_COUNT > 0) contourAmortization(contour);
    centerContour(contour);

    std::vector<Vector2> vec_pixels;
    vec_pixels.reserve(contour.size());
    for (const auto& pix : contour){
        vec_pixels.push_back(Vector2(static_cast<float>(pix.real()) + offset.x, static_cast<float>(pix.imag()) + offset.y));
    }

    auto epicycles = computeEpicycles(contour);
    
    // --- ПРЕАЛЛОКАЦИЯ БУФЕРОВ ДЛЯ ГОРЯЧЕГО ЦИКЛА ---
    VectorArray frame_vectors;
    JointArray frame_joints;

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
    bool showOutline = true;
    double speedMult = 0.2; // Начальная скорость

    while (!WindowShouldClose()) {
        // Управление скоростью стрелочками ВВЕРХ / ВНИЗ
        if (IsKeyPressed(KEY_UP)) speedMult += 0.02;
        if (IsKeyPressed(KEY_DOWN) && speedMult > 0.03) speedMult -= 0.02;
        if (IsKeyPressed(KEY_S)) speedMult = 0.004;

        // Физика
        calculateFrame(epicycles, t, frame_vectors, frame_joints);
        if (showTrail && (speedMult * trail.size() < pixels_count)) {
            trail.push_back({ static_cast<float>(frame_joints.back().pos.real()) + offset.x, static_cast<float>(frame_joints.back().pos.imag()) + offset.y });
        }

        t += dt * speedMult;
        t -= (int)t;

        // Отрисовка
        BeginDrawing();
        ClearBackground({20, 20, 20, 255}); // Темный фон

        // Рисуем эпициклы
        for (size_t i = 0; i < NUM_EPICYCLES; ++i) {
            Vector2 p1 = { static_cast<float>(frame_joints[i].pos.real()) + offset.x, static_cast<float>(frame_joints[i].pos.imag()) + offset.y };
            Vector2 p2 = { static_cast<float>(frame_joints[i+1].pos.real()) + offset.x, static_cast<float>(frame_joints[i+1].pos.imag()) + offset.y };
            float radius = static_cast<float>(frame_joints[i].radius);

            // Рисуем круги только если радиус > 0.5 пикселя, иначе их всё равно не видно
            if (showCircles && radius > 0.5f) {
                DrawCircleLines(p1.x, p1.y, radius, Fade(WHITE, 0.15f));
            }
            if (showLines) {
                DrawLineV(p1, p2, WHITE);
            }
        }

        // Рисуем след (polyline)
        if (showTrail && trail.size() > 1) {
            DrawLineStrip(trail.data(), trail.size(), YELLOW);
        }

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