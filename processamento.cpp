#include "PROCESSAMENTO DE IMAGEM/image.hpp"
#include "PROCESSAMENTO DE IMAGEM/edgeDetect.hpp"

int main() {
    Image img = Image::loadGray("../PROCESSAMENTO DE IMAGEM/images/neymar.jpg");
    Image bordas_sobel = alg1_sobel(img, 100);
    Image bordas_log   = alg2_log(img, 20);

    bordas_sobel.save("../PROCESSAMENTO DE IMAGEM/saida_sobel.jpg");
    bordas_log.save("../PROCESSAMENTO DE IMAGEM/saida_log.jpg");
    return 0;
}