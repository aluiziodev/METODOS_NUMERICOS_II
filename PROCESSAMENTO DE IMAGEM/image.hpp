#ifndef IMAGE_HPP
#define IMAGE_HPP
#define STB_IMAGE_IMPLEMENTATION
#include "../libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../libs/stb_image_write.h"
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>



struct Image {
    int width, height, channels;
    std::vector<unsigned char> data; // grayscale: 1 canal

    static Image loadGray(const std::string &path) {
        Image img;
        unsigned char* raw = stbi_load(path.c_str(), &img.width, &img.height, &img.channels, 1);
        if (!raw) {
            std::cerr << "ERRO ao carregar imagem: " << path << std::endl;
            std::cerr << "Motivo: " << stbi_failure_reason() << std::endl;
            std::exit(1);
        }
        img.channels = 1;
        img.data.assign(raw, raw + img.width * img.height);
        stbi_image_free(raw);
        return img;
    }

    unsigned char get(int i, int j) const {
        // clamp nas bordas (extensão)
        int ii = std::min(std::max(i, 0), height - 1);
        int jj = std::min(std::max(j, 0), width - 1);
        return data[ii * width + jj];
    }

    void save(const std::string &path) const {
        int ok = stbi_write_jpg(path.c_str(), width, height, 1, data.data(), 90);
        if (!ok) {
            std::cerr << "ERRO ao salvar imagem: " << path << std::endl;
        } else {
            std::cout << "Imagem salva em: " << path << std::endl;
        }
    }
};
#endif