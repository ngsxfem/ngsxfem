#
## Cube and Spheres
#
algebraic3d

# a cube
solid cube = plane (-1,-1,-1; 0, 0,-1)
         and plane (-1,-1,-1; 0,-1, 0)
         and plane (-1,-1,-1;-1, 0, 0)
         and plane ( 1, 1, 1; 0, 0, 1)
         and plane ( 1, 1, 1; 0, 1, 0)
         and plane ( 1, 1, 1; 1, 0, 0);

tlo cube;
