#ifndef CONVOLVE_HPP
#define CONVOLVE_HPP
#include "image.hpp"
#include "edgeKernels.hpp"

Image convolve(const Image &img, const Kernel3x3 &k, double norm = 1.0) {
    Image out = img;
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {
            double sum = 0.0;
            for (int di = -1; di <= 1; di++)
                for (int dj = -1; dj <= 1; dj++)
                    sum += k[di+1][dj+1] * img.get(i+di, j+dj);
            sum /= norm;
            out.data[i*img.width + j] = (unsigned char) std::min(255.0, std::max(0.0, sum));
        }
    }
    return out;
}
#endif