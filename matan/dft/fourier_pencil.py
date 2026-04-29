import pygame
import sys

# Настройки
IMAGE_PATH = "photo.png" # Положи фото в папку со скриптом
OUTPUT_FILE = "points.txt"

pygame.init()
img = pygame.image.load(IMAGE_PATH)
screen = pygame.display.set_mode(img.get_size())
pygame.display.set_caption("Обведи контур (не отрывая мышку) и закрой окно")

points = []
drawing = False

clock = pygame.time.Clock()

while True:
    screen.blit(img, (0, 0))
    
    # Рисуем уже набранные точки, чтобы видеть прогресс
    if len(points) > 1:
        pygame.draw.lines(screen, (255, 0, 0), False, points, 2)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            # Сохраняем и выходим
            with open(OUTPUT_FILE, "w") as f:
                for p in points:
                    f.write(f"{p[0]} {p[1]}\n")
            pygame.quit()
            sys.exit()
            
        if event.type == pygame.MOUSEBUTTONDOWN:
            drawing = True
            points = [] # Начинаем новый контур
        
        if event.type == pygame.MOUSEBUTTONUP:
            drawing = False
            
        if event.type == pygame.MOUSEMOTION and drawing:
            pos = pygame.mouse.get_pos()
            # Добавляем точку, только если мышка сдвинулась (чтобы не дублировать)
            if not points or (pygame.math.Vector2(pos) - pygame.math.Vector2(points[-1])).length() > 2:
                points.append(pos)

    pygame.display.flip()
    clock.tick(60)