#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <numeric>
#include <print>
#include <ranges>
#include <string>
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

void calculateFrame(const std::vector<Epicycle>& epicycles, double t, std::vector<Complex>& vectors, std::vector<Joint>& joints) noexcept{
    if (epicycles.empty()) return;
    // [[assume(epicycles.size() == vectors.size())]];
    // [[assume(epicycles.size() + 1 == joints.size())]];

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
    constexpr int M = 100;
    constexpr int width = 1280;
    constexpr int height = 720;
    constexpr Vector2 offset = { width / 2.0f, height / 2.0f };

    const std::vector<std::pair<int, int>> raw_pixels = {{287, 714}, {279, 704}, {264, 682}, {243, 648}, {217, 608}, {194, 565}, {181, 528}, {181, 495}, {191, 465}, {213, 438}, {251, 414}, {300, 396}, {352, 385}, {404, 384}, {450, 394}, {478, 410}, {496, 436}, {513, 473}, {527, 503}, {542, 530}, {555, 548}, {561, 556}, {565, 557}, {573, 555}, {585, 550}, {603, 543}, {623, 533}, {644, 521}, {664, 507}, {687, 493}, {710, 480}, {734, 468}, {753, 459}, {771, 456}, {788, 456}, {804, 460}, {821, 470}, {836, 483}, {850, 497}, {862, 514}, {869, 530}, {875, 548}, {877, 566}, {876, 580}, {873, 598}, {868, 613}, {859, 628}, {847, 642}, {830, 659}, {810, 674}, {789, 689}, {770, 703}, {750, 716}, {729, 725}, {709, 732}, {690, 736}, {670, 736}, {652, 735}, {637, 732}, {623, 725}, {611, 717}, {601, 710}, {593, 703}, {585, 697}, {577, 691}, {570, 687}, {564, 684}, {556, 683}, {548, 683}, {539, 686}, {528, 691}, {516, 699}, {505, 710}, {497, 722}, {489, 738}, {483, 754}, {478, 772}, {477, 790}, {479, 808}, {481, 824}, {484, 841}, {487, 855}, {491, 870}, {494, 882}, {497, 894}, {498, 901}, {498, 906}, {495, 909}, {490, 910}, {482, 910}, {474, 911}, {466, 910}, {457, 908}, {448, 905}, {439, 902}, {431, 900}, {421, 897}, {413, 895}, {404, 893}, {395, 893}, {385, 893}, {377, 892}, {367, 890}, {357, 887}, {345, 885}, {335, 885}, {323, 885}, {314, 888}, {305, 891}, {297, 895}, {289, 901}, {283, 909}, {275, 920}, {266, 929}, {258, 936}, {252, 940}, {248, 943}, {244, 945}, {241, 947}, {237, 950}, {230, 952}, {223, 954}, {214, 954}, {203, 954}, {191, 953}, {177, 951}, {162, 947}, {148, 944}, {139, 939}, {135, 935}, {132, 930}, {128, 924}, {124, 919}, {120, 913}, {116, 908}, {113, 902}, {111, 896}, {107, 890}, {104, 883}, {102, 877}, {101, 870}, {102, 860}, {105, 851}, {108, 840}, {113, 829}, {118, 820}, {123, 813}, {128, 807}, {133, 803}, {137, 800}, {141, 797}, {145, 795}, {151, 794}, {160, 794}, {171, 794}, {185, 794}, {197, 794}, {208, 795}, {218, 796}, {229, 798}, {236, 799}, {244, 799}, {255, 799}, {267, 799}, {280, 799}, {291, 799}, {303, 798}, {314, 797}, {324, 795}, {332, 793}, {337, 790}, {343, 785}, {350, 778}, {358, 769}, {364, 757}, {370, 745}, {376, 732}, {381, 721}, {385, 710}, {387, 699}, {389, 688}, {391, 676}, {393, 665}, {394, 653}, {395, 642}, {396, 631}, {397, 619}, {398, 606}, {399, 597}, {399, 588}, {399, 579}, {399, 569}, {400, 559}, {400, 548}, {399, 538}, {396, 530}, {392, 524}, {388, 521}, {385, 518}, {382, 516}, {377, 514}, {372, 512}, {369, 511}, {367, 509}, {365, 507}, {361, 505}, {357, 503}, {350, 500}, {344, 498}, {340, 497}, {336, 497}, {333, 499}, {332, 502}, {330, 507}, {330, 512}, {330, 517}, {332, 523}, {335, 532}, {337, 541}, {338, 550}, {340, 558}, {342, 566}, {343, 572}, {343, 578}, {343, 583}, {343, 587}, {341, 588}, {338, 588}, {334, 587}, {330, 586}, {326, 585}, {322, 583}, {315, 578}, {309, 572}, {302, 564}, {298, 556}, {294, 548}, {289, 538}, {285, 530}, {281, 524}, {280, 520}, {282, 515}, {286, 509}, {292, 502}, {299, 496}, {306, 489}, {314, 483}, {322, 478}, {331, 472}, {343, 465}, {353, 459}, {363, 455}, {371, 452}, {374, 451}, {377, 452}, {381, 456}, {387, 463}, {396, 472}, {407, 484}, {419, 497}, {429, 508}, {437, 518}, {446, 525}, {451, 530}, {455, 534}, {458, 537}, {462, 541}, {468, 547}, {477, 555}, {486, 564}, {495, 573}, {503, 580}, {509, 587}, {518, 595}, {528, 604}, {540, 615}, {552, 625}, {562, 631}, {572, 637}, {581, 641}, {591, 644}, {604, 648}, {618, 651}, {631, 653}, {644, 654}, {657, 654}, {669, 653}, {682, 653}, {691, 653}, {699, 653}, {708, 655}, {719, 657}, {731, 657}, {743, 657}, {752, 656}, {761, 654}, {767, 652}, {772, 649}, {775, 646}, {778, 642}, {779, 638}, {780, 633}, {780, 626}, {780, 617}, {780, 608}, {780, 602}, {778, 597}, {775, 594}, {772, 594}, {770, 593}, {764, 593}, {761, 594}, {756, 594}, {750, 594}, {741, 595}, {732, 595}, {724, 595}, {715, 595}, {706, 594}, {699, 594}, {689, 592}, {679, 591}, {669, 590}, {658, 588}, {648, 586}, {638, 585}, {629, 585}, {620, 585}, {613, 585}, {606, 585}, {598, 585}, {591, 584}, {585, 584}, {580, 583}, {574, 583}, {567, 582}, {560, 581}, {553, 579}, {547, 579}, {543, 577}, {538, 575}, {533, 572}, {528, 569}, {524, 567}, {521, 565}, {519, 563}, {515, 559}, {509, 553}, {503, 547}, {497, 541}, {492, 536}, {487, 531}, {480, 525}, {475, 519}, {471, 514}, {467, 508}, {463, 503}, {460, 498}, {456, 494}, {453, 489}, {449, 485}, {445, 480}, {441, 476}, {437, 472}, {434, 469}, {430, 466}, {427, 464}, {423, 460}, {418, 457}, {413, 454}, {407, 451}, {403, 448}, {399, 446}, {394, 444}, {390, 441}, {384, 438}, {376, 435}, {368, 432}, {361, 431}, {356, 431}, {351, 431}, {346, 433}, {343, 434}, {339, 435}, {333, 437}, {328, 438}, {324, 440}, {319, 442}, {316, 443}, {312, 446}, {307, 450}, {302, 454}, {298, 458}, {294, 464}, {291, 469}, {288, 475}, {285, 481}, {282, 486}, {280, 492}, {279, 498}, {278, 503}, {278, 509}, {277, 514}, {277, 520}, {276, 525}, {275, 532}, {275, 540}, {274, 547}, {273, 553}, {273, 558}, {273, 561}, {273, 565}, {273, 569}, {273, 573}, {273, 578}, {272, 584}, {272, 589}, {273, 594}, {274, 599}, {276, 605}, {277, 610}, {278, 614}, {280, 618}, {282, 621}, {284, 625}, {287, 629}, {290, 634}, {293, 639}, {296, 642}, {298, 644}, {300, 647}, {302, 650}, {304, 654}, {307, 658}, {310, 662}, {312, 665}, {314, 667}, {316, 669}, {318, 673}, {321, 676}, {324, 679}, {325, 681}, {326, 683}, {327, 685}, {328, 689}, {329, 693}, {330, 696}, {331, 699}, {332, 701}, {333, 705}, {334, 708}, {334, 711}, {334, 714}, {334, 717}, {333, 720}, {333, 723}, {332, 726}, {330, 730}, {328, 734}, {325, 739}, {322, 743}, {319, 747}, {317, 750}, {314, 752}, {310, 755}, {304, 759}, {298, 762}, {291, 766}, {284, 770}, {276, 773}, {269, 775}, {260, 777}, {252, 778}, {244, 780}, {234, 780}, {228, 781}, {221, 781}, {214, 781}, {207, 781}, {201, 780}, {194, 779}, {188, 778}, {183, 778}, {179, 778}, {175, 777}, {170, 776}, {165, 773}, {161, 771}, {157, 769}, {153, 765}, {148, 761}, {143, 757}, {138, 753}, {134, 749}, {131, 745}, {126, 741}, {123, 737}, {120, 732}, {118, 729}, {115, 726}, {113, 722}, {110, 717}, {107, 712}, {104, 707}, {101, 702}, {100, 699}, {98, 695}, {96, 690}, {95, 686}, {95, 682}, {95, 679}, {96, 677}, {97, 675}, {99, 673}, {101, 672}, {104, 672}, {109, 673}, {115, 675}, {120, 677}, {127, 680}, {135, 684}, {142, 688}, {149, 692}, {154, 695}, {159, 698}, {163, 701}, {168, 703}, {173, 705}, {178, 707}, {182, 709}, {186, 711}, {191, 714}, {197, 717}, {204, 721}, {210, 724}, {215, 726}, {219, 728}, {222, 729}, {226, 730}, {230, 731}, {233, 733}, {237, 733}, {240, 733}, {243, 733}, {247, 733}, {251, 733}, {254, 733}, {257, 733}, {261, 733}, {265, 733}, {268, 733}, {271, 732}, {273, 731}, {275, 730}, {278, 728}, {280, 726}, {281, 724}, {282, 722}, {284, 721}, {285, 719}, {286, 717}, {287, 715}, {288, 713}};
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
    bool showOutline = true;
    double speedMult = 0.1; // Начальная скорость

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
        DrawLineStrip(vec_pixels.data(), vec_pixels.size(), Fade(WHITE, 0.05f));

        // --- ОТРИСОВКА UI ---
        DrawToggleButton({20, 20, 140, 40}, "Circles", showCircles);
        DrawToggleButton({170, 20, 140, 40}, "Lines", showLines);
        DrawToggleButton({320, 20, 140, 40}, "Curve", showTrail);
        DrawToggleButton({470, 20, 140, 40}, "Curve", showOutline);
        DrawClearButton({620, 20, 140, 40}, "Clear", trail);
        
        // Текст со скоростью
        DrawText(TextFormat("Speed: %.2fx (Arrows UP/DOWN)", speedMult), 780, 30, 20, LIGHTGRAY);
        DrawText(TextFormat("Points: %zu | M: %d", trail.size(), M), 20, height - 30, 20, DARKGRAY);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}