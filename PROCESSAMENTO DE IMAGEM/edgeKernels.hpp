#ifndef EDGE_KERNELS_HPP
#define EDGE_KERNELS_HPP
#include <array>

using Kernel3x3 = std::array<std::array<double,3>,3>;

// Sobel Gx e Gy (equivalem à 1ª derivada central O(h) aplicada em 2D)
Kernel3x3 sobelX() {
    return {{ {-1,0,1}, {-2,0,2}, {-1,0,1} }};
}
Kernel3x3 sobelY() {
    return {{ {-1,-2,-1}, {0,0,0}, {1,2,1} }};
}

// Laplaciano (soma das 2ª derivadas centrais em x e y)
Kernel3x3 laplace() {
    return {{ {0,1,0}, {1,-4,1}, {0,1,0} }};
}

// Kernel gaussiano 3x3, sigma=1.0 (para suavização antes da detecção)
Kernel3x3 gauss3x3() {
    return {{ {1,2,1}, {2,4,2}, {1,2,1} }}; // depois normaliza dividindo por 16
}
#endif