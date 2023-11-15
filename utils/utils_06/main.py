# Pygame
import pygame as pg
import sys
sys.path.insert(0, '')

# local classes
from utils.utils_06.boid import Boid

WIDTH = 800
HEIGHT = 800

def run(WIDTH, HEIGHT, num_boids, a, c, s):
    pg.init()
    screen = pg.display.set_mode((WIDTH, HEIGHT))
    screen.fill("#FFFAFA")
    pg.display.set_caption("Boids")
    clock = pg.time.Clock()

    # Create boids, hoiks and obstacles
    boids = [Boid(WIDTH, HEIGHT) for i in range(num_boids)]

    running = True
    while running:
        for event in pg.event.get():
            if event.type == pg.QUIT:
                running = False
        screen.fill((255,250,250))

        # Boid loop 
        for boid in boids: 
            boid.draw(screen)
            boid.update(boids, a, c, s)

        # Update the screen
        pg.display.flip()
        clock.tick(60)
    pg.quit()

if __name__ == "__main__":
    num_boids = int(sys.argv[1])
    a = float(sys.argv[2])
    c = float(sys.argv[3])
    s = float(sys.argv[4])

    run(WIDTH, HEIGHT, num_boids, a, c, s)
