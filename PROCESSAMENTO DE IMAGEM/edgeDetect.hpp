#ifndef EDGE_DETECT_HPP
#define EDGE_DETECT_HPP
#include "image.hpp"
#include "convolve.hpp"
#include <cmath>

Image alg1_sobel(const Image &img, unsigned char threshold) {
    Image blurred = convolve(img, gauss3x3(), 16.0);
    Image gx = convolve(blurred, sobelX());
    Image gy = convolve(blurred, sobelY());

    Image out = img;
    for (size_t p = 0; p < img.data.size(); p++) {
        double mag = std::sqrt(gx.data[p]*gx.data[p] + gy.data[p]*gy.data[p]);
        out.data[p] = (mag > threshold) ? 255 : 0;
    }
    return out;
}

Image alg2_log(const Image &img, unsigned char tolerance) {
    Image blurred = convolve(img, gauss3x3(), 16.0);
    Image log_img = convolve(blurred, laplace());

    Image out = img;
    for (size_t p = 0; p < img.data.size(); p++) {
        out.data[p] = (log_img.data[p] > tolerance) ? 255 : 0;
    }
    return out;
}
#endif