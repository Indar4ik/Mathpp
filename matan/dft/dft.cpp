#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <numeric>
#include <print>
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

    double sum_re = 0.0;
    double sum_im = 0.0;

    for (const auto& pt : raw_points){
        double re = static_cast<double>(pt.first);
        double im = static_cast<double>(pt.second);
        contour.emplace_back(re, im);
        sum_re += re;
        sum_im += im;
    }

    Complex center(sum_re / n, sum_im / n);

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

    #pragma clang unroll(enable)
    for (int n = -M; n <= M; ++n) {
        double sum_re = 0.0;
        double sum_im = 0.0;

        #pragma clang loop vectorize(enable)
        for (size_t k = 0; k < N; ++k) {
            double angle = pi2_N * n * static_cast<double>(k);
            
            double c = std::cos(angle);
            double s = std::sin(angle);
            
            double re = contour[k].real();
            double im = contour[k].imag();

            sum_re += re * c - im * s;
            sum_im += re * s + im * c;
        }
        
        Complex c_n(sum_re / N, sum_im / N);
        // std::println("c_{}: {}", n, std::abs(c_n));
        epicycles.push_back({static_cast<double>(n), c_n, std::abs(c_n)});
    }

    std::ranges::sort(epicycles, std::ranges::greater{},[](const Epicycle& e) { return std::norm(e.coef); });

    return epicycles;
}

void calculateFrame(const std::vector<Epicycle>& epicycles, double t, std::vector<Complex>& vectors, std::vector<Joint> joints) noexcept{
    if (epicycles.empty()) return;
    [[assume(epicycles.size() == vectors.size())]];
    [[assume(epicycles.size() + 1 == joints.size())]];

    const size_t n = epicycles.size();

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

    Complex cur_pos = {0.0, 0.0};
    for (size_t i = 0; i < n; ++i){
        joints[i] = {cur_pos, epicycles[i].radius};
        cur_pos += vectors[i];
    }
    joints[n] = {cur_pos, 0.0};
}

std::vector<std::pair<int, int>> generateHeart() {
    std::vector<std::pair<int, int>> points;
    for (double t = 0; t < 2 * std::numbers::pi; t += 0.02) {
        // Формула сердечка
        double x = 240 * std::pow(std::sin(t), 3);
        double y = 15 * (13 * std::cos(t) - 5 * std::cos(2*t) - 2 * std::cos(3*t) - std::cos(4*t));
        // Масштабируем
        std::print("({:.4f}, {:.4f}) ", x, y);
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

bool DrawClearButton(Rectangle bounds, const char* text, std::vector<Vector2>& vec) {
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
    constexpr int M = 50;
    constexpr int width = 1280;
    constexpr int height = 720;
    constexpr Vector2 offset = { width / 2.0f, height / 2.0f };

    const std::vector<std::pair<int, int>> raw_pixels = {{330, 413}, {328, 414}, {326, 415}, {323, 416}, {321, 417}, {319, 418}, {317, 421}, {316, 427}, {316, 434}, {316, 442}, {316, 449}, {317, 457}, {321, 466}, {326, 475}, {332, 483}, {338, 489}, {344, 493}, {351, 496}, {360, 499}, {375, 501}, {387, 503}, {397, 504}, {405, 505}, {411, 505}, {417, 505}, {424, 502}, {433, 497}, {441, 492}, {447, 487}, {451, 481}, {454, 474}, {457, 467}, {460, 460}, {461, 454}, {461, 446}, {460, 440}, {457, 433}, {453, 425}, {449, 418}, {444, 412}, {438, 407}, {432, 403}, {426, 400}, {419, 398}, {412, 395}, {405, 393}, {400, 391}, {395, 391}, {390, 391}, {385, 391}, {379, 393}, {373, 395}, {368, 397}, {365, 399}, {362, 400}, {359, 401}, {357, 402}, {354, 403}, {351, 404}, {349, 405}, {346, 405}, {344, 406}, {340, 406}, {337, 407}, {335, 408}, {331, 409}, {329, 410}, {334, 410}, {336, 409}, {339, 408}, {342, 408}, {345, 407}, {348, 406}, {352, 406}, {355, 405}, {359, 404}, {363, 403}, {367, 402}, {371, 402}, {374, 401}, {380, 401}, {385, 401}, {388, 401}, {391, 401}, {394, 400}, {399, 400}, {403, 400}, {406, 399}, {410, 400}, {413, 400}, {416, 400}, {420, 401}, {423, 402}, {425, 404}, {428, 405}, {430, 407}, {434, 409}, {437, 412}, {440, 415}, {442, 418}, {444, 420}, {445, 422}, {445, 425}, {447, 428}, {448, 431}, {449, 436}, {450, 438}, {450, 441}, {451, 443}, {452, 446}, {453, 448}, {454, 450}, {455, 452}, {456, 454}, {457, 456}, {461, 456}, {466, 458}, {470, 459}, {473, 460}, {477, 460}, {480, 460}, {483, 460}, {486, 460}, {490, 460}, {494, 461}, {498, 462}, {501, 463}, {503, 464}, {505, 465}, {508, 467}, {511, 468}, {514, 470}, {516, 471}, {518, 472}, {521, 474}, {524, 475}, {526, 476}, {528, 477}, {530, 478}, {532, 479}, {536, 480}, {538, 479}, {540, 478}, {541, 476}, {543, 473}, {545, 471}, {548, 468}, {550, 466}, {552, 465}, {555, 463}, {559, 461}, {562, 459}, {567, 457}, {572, 454}, {576, 452}, {580, 451}, {584, 450}, {588, 449}, {592, 449}, {596, 449}, {600, 449}, {605, 449}, {610, 449}, {615, 450}, {620, 451}, {623, 452}, {628, 454}, {633, 456}, {638, 459}, {644, 461}, {649, 463}, {655, 466}, {660, 468}, {664, 472}, {670, 475}, {674, 478}, {678, 481}, {681, 484}, {683, 487}, {686, 491}, {688, 495}, {690, 499}, {691, 503}, {692, 507}, {693, 512}, {693, 515}, {693, 520}, {691, 523}, {689, 526}, {685, 530}, {680, 534}, {675, 537}, {670, 540}, {666, 543}, {661, 545}, {657, 547}, {653, 548}, {649, 550}, {644, 551}, {640, 553}, {635, 554}, {630, 556}, {626, 556}, {622, 556}, {617, 556}, {610, 556}, {604, 555}, {598, 555}, {591, 553}, {586, 551}, {582, 549}, {576, 547}, {571, 543}, {565, 539}, {561, 536}, {557, 533}, {554, 530}, {551, 527}, {548, 524}, {546, 521}, {544, 519}, {542, 516}, {541, 513}, {539, 511}, {538, 508}, {537, 505}, {537, 501}, {536, 498}, {536, 494}, {536, 490}, {537, 488}, {538, 486}, {540, 485}, {540, 482}, {541, 480}, {540, 477}};
    const double pixels_count = static_cast<double>(raw_pixels.size());

    auto contour = prepareContour(raw_pixels);

    std::vector<Vector2> vec_pixels;
    vec_pixels.reserve(raw_pixels.size());
    for (const auto& pix : contour){
        vec_pixels.push_back(Vector2(static_cast<float>(pix.real()) + offset.x, static_cast<float>(pix.imag()) + offset.y));
    }

    auto epicycles = computeEpicycles(contour, M);
    
    // --- ПРЕАЛЛОКАЦИЯ БУФЕРОВ ДЛЯ ГОРЯЧЕГО ЦИКЛА ---
    std::vector<Complex> frame_vectors(epicycles.size());
    std::vector<Joint> frame_joints(epicycles.size() + 1);

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
        if (IsKeyPressed(KEY_DOWN) && speedMult > 0.06) speedMult -= 0.05;

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
        for (size_t i = 0; i < frame_joints.size() - 1; ++i) {
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
        DrawLineStrip(vec_pixels.data(), vec_pixels.size(), Fade(WHITE, 0.1f));

        // --- ОТРИСОВКА UI ---
        DrawToggleButton({20, 20, 140, 40}, "Circles", showCircles);
        DrawToggleButton({170, 20, 140, 40}, "Lines", showLines);
        DrawToggleButton({320, 20, 140, 40}, "Curve", showTrail);
        DrawClearButton({470, 20, 140, 40}, "Clear", trail);
        
        // Текст со скоростью
        DrawText(TextFormat("Speed: %.2fx (Arrows UP/DOWN)", speedMult), 630, 30, 20, LIGHTGRAY);
        DrawText(TextFormat("Points: %zu | M: %d", trail.size(), M), 20, height - 30, 20, DARKGRAY);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}