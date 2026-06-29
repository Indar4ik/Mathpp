import os
import sys
import math

# --- Wayland fix ---
# Под Wayland SDL2 отдаёт координаты мыши в масштабированном пространстве, из-за чего
# штрих перекашивается ("криво"). Форсируем XWayland (как на Windows) ДО импорта pygame.
# Если у вас чистый X11 — строка ничего не ломает. Можно отключить, выставив PENCIL_NATIVE=1.
if sys.platform.startswith("linux") and not os.environ.get("PENCIL_NATIVE"):
    os.environ.setdefault("SDL_VIDEODRIVER", "x11")

import pygame

# --- Настройки ---
IMAGE_PATH = "photo.png"   # Фон для обводки (можно отсутствующий -> пустой холст)
OUTPUT_FILE = "points.txt"
STEP = 3.0                 # Желаемый шаг между точками в пикселях (плотность контура)
FPS = 240

pygame.init()

try:
    img = pygame.image.load(IMAGE_PATH)   # .convert() — только ПОСЛЕ set_mode
    size = img.get_size()
except (pygame.error, FileNotFoundError):
    img = None
    size = (1280, 720)

screen = pygame.display.set_mode(size)
if img is not None:
    img = img.convert()                   # ускоряет блит, требует готового видеорежима
pygame.display.set_caption("ЛКМ — рисовать | C — очистить | S/Enter — сохранить | Esc — выход")
font = pygame.font.SysFont("monospace", 18)
clock = pygame.time.Clock()

# Полупрозрачный слой затемнения фона: чем темнее фото, тем заметнее контур.
# Регулируется стрелками ВВЕРХ/ВНИЗ.
dim = 120
dim_overlay = pygame.Surface(size, pygame.SRCALPHA)

points = []            # список (float, float)
drawing = False
saved_flash = 0        # кадры, в течение которых показываем "Saved"


def hud(text, color):
    # convert_alpha() обязателен: блит текста без него на convert()-фон в формате дисплея
    # x11 даёт артефакт "сплошной прямоугольник". Приводим к формату дисплея с альфой.
    return font.render(text, True, color).convert_alpha()


def add_sampled(pos):
    """Добавляет точку pos, доинтерполируя промежуточные с шагом ~STEP,
    чтобы быстрые движения не оставляли разрывов (платформо-независимо)."""
    if not points:
        points.append((float(pos[0]), float(pos[1])))
        return
    lx, ly = points[-1]
    dx, dy = pos[0] - lx, pos[1] - ly
    dist = math.hypot(dx, dy)
    if dist < STEP:
        return
    n = int(dist / STEP)
    for i in range(1, n + 1):
        f = (i * STEP) / dist
        points.append((lx + dx * f, ly + dy * f))


def save():
    with open(OUTPUT_FILE, "w") as f:
        for x, y in points:
            f.write(f"{int(round(x))} {int(round(y))}\n")


while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            if points:
                save()
            pygame.quit()
            sys.exit()

        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_ESCAPE:
                pygame.quit()
                sys.exit()
            if event.key in (pygame.K_s, pygame.K_RETURN, pygame.K_KP_ENTER):
                if points:
                    save()
                    saved_flash = FPS  # ~1 секунда
            if event.key == pygame.K_c:
                points = []
            if event.key == pygame.K_UP:
                dim = min(dim + 20, 240)
            if event.key == pygame.K_DOWN:
                dim = max(dim - 20, 0)

        if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
            drawing = True
            points = []                       # новый контур одним штрихом
            add_sampled(event.pos)

        if event.type == pygame.MOUSEBUTTONUP and event.button == 1:
            drawing = False

    # Семплируем позицию КАЖДЫЙ кадр, а не по MOUSEMOTION — равномерно и без рваных пропусков
    if drawing:
        add_sampled(pygame.mouse.get_pos())

    # --- Отрисовка ---
    if img is not None:
        screen.blit(img, (0, 0))
        if dim > 0:                       # затемняем фон, чтобы контур был заметнее
            dim_overlay.fill((0, 0, 0, dim))
            screen.blit(dim_overlay, (0, 0))
    else:
        screen.fill((20, 20, 20))

    if len(points) > 1:
        # aalines — сглаженная линия, выглядит ровнее ломаной
        pygame.draw.aalines(screen, (255, 60, 60), False, points)
        # маркер старта и пунктир замыкания (контур в C++ замыкается автоматически)
        sx, sy = points[0]
        pygame.draw.circle(screen, (60, 200, 60), (int(sx), int(sy)), 5, 2)
        pygame.draw.line(screen, (80, 80, 80), points[0], points[-1], 1)

    # HUD на полупрозрачной подложке — читается над любым фоном
    info = f"points: {len(points)}   dim: {dim}   [C]lear [S]ave [^/v]dim [Esc]quit"
    txt = hud(info, (255, 255, 0))
    panel = pygame.Surface((txt.get_width() + 12, txt.get_height() + 8), pygame.SRCALPHA)
    panel.fill((0, 0, 0, 140))
    screen.blit(panel, (6, 6))
    screen.blit(txt, (12, 10))
    if saved_flash > 0:
        screen.blit(hud(f"Saved -> {OUTPUT_FILE}", (60, 220, 60)), (12, 34))
        saved_flash -= 1

    pygame.display.flip()
    clock.tick(FPS)
