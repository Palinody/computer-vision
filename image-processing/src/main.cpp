#include "headers/PngParser.h"
#include "headers/ImageProcessing.h"
#include "headers/Timer.h"
#include "headers/PRNG.h"
#include "headers/BitMap.h"

#include "headers/MipMap.h"

#include <stdio.h>
#include <stdlib.h> //itoa atoi
#include <cmath> // gamma correction power
#include <algorithm> // std::max(a, b)
/**
 * Compilation: 
 *      g++-10 -std=c++20 -Wall -O3 -ffast-math -march=native -o main main.cpp -lpng
 * Tranformations:
 *      github avg simple method: https://stackoverflow.com/questions/9570895/image-downscaling-algorithm
 *      Lanczos resampling: https://en.wikipedia.org/wiki/Lanczos_resampling#Advantages
 *      Bicubic interpolation: https://en.wikipedia.org/wiki/Bicubic_interpolation
 *      Bilinear interpolation: https://en.wikipedia.org/wiki/Bilinear_interpolation
 *      Spline interpolation: https://en.wikipedia.org/wiki/Spline_interpolation
 * 
 * libpng source code: https://github.com/LuaDist/libpng
 * 
 * gif maker one header code https://github.com/charlietangora/gif-h
 * wiki: https://en.wikipedia.org/wiki/GIF#Example_GIF_file
 * 
 * Animated Portable Network Graphics (APNG): https://en.wikipedia.org/wiki/APNG
 * code: https://sourceforge.net/p/libpng-apng/code/ci/master/tree/
 * 
 * libpng manual: http://www.libpng.org/pub/png/libpng-1.2.5-manual.html
*/

// init image_manager
// get image_manager dims info and init pixelmap
// pass pixelmap to 

void downscale_ex(){
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/processed.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);


    PixelMap<png_byte> pixelmap_blurred = pixelmap;
    // orig / 2
    PixelMap<png_byte> pixelmap2(pixelmap.getHeight()/2, pixelmap.getWidth()/2, pixelmap.getChannels());
    
    PixelMap<png_byte> pixelmap2_blurred = pixelmap2;
    // orig / 4
    PixelMap<png_byte> pixelmap4(pixelmap2.getHeight()/2, pixelmap2.getWidth()/2, pixelmap2.getChannels());
    
    PixelMap<png_byte> pixelmap4_blurred = pixelmap4;
    // orig / 8
    PixelMap<png_byte> pixelmap8(pixelmap4.getHeight()/2, pixelmap4.getWidth()/2, pixelmap4.getChannels());

    double sd = 0.2;
    size_t k_w = std::ceil(sd * 6);
    Kernel1D<float> ker(k_w, dims.channels, sd);


    fft(pixelmap_blurred, pixelmap, ker);
    rescale::box(pixelmap2, pixelmap_blurred);

    fft(pixelmap2_blurred, pixelmap2, ker);
    rescale::box(pixelmap4, pixelmap2_blurred);
    
    fft(pixelmap4_blurred, pixelmap4, ker);
    rescale::box(pixelmap8, pixelmap4_blurred);
    
    write_png_file(to_file, pixelmap8);
}

void blurr_ex(){
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/blurred.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    PixelMap<png_byte> inter(src.getHeight(), src.getWidth(), src.getChannels());
    
    read_png_file(from_file, src);

    double sd = 1.1;
    size_t k_w = std::ceil(sd * 6);
    Kernel1D<float> ker(k_w, src.getChannels(), sd);

    /*
    blurr(inter, src, ker);
    inter = inter.transpose();
    PixelMap<png_byte> output(inter.getHeight(), inter.getWidth(), inter.getChannels());
    blurr(output, inter, ker);
    output = output.transpose();
    */
    PixelMap<png_byte> output(dims.height, dims.width, dims.channels);
    fft(output, src, ker);

    printf("channels: %ld\n", output.getChannels());
    write_png_file(to_file, output);
}

float nearestNeighbor(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/nearestNeighbor.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);

    PixelMap<png_byte> dest(
        dims.height*10,
        dims.width*10,
        dims.channels
    );
    timer.reset();
    rescale::nearestNeighbor(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (nearest neighbor): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float nearestNeighbor_bitmap_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/nearestNeighbor_bitmap_rgba.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGBA src(pixelmap);

    BitMapRGBA dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::nearestNeighbor(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (nearest neighbor - bitmap - rgba): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float nearestNeighbor_bitmap_rgb(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/nearestNeighbor_bitmap_rgb.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGB src(pixelmap);

    BitMapRGB dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::nearestNeighbor(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (nearest neighbor - bitmap - rgb): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bilinear(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bilinear.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);

    PixelMap<png_byte> dest(
        dims.height*10,
        dims.width*10,
        dims.channels
    );
    timer.reset();
    rescale::bilinear(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bilinear): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bilinear_bitmap_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bilinear_bitmap_rgba.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGBA src(pixelmap);

    BitMapRGBA dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::bilinear(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bilinear - bitmap - rgba): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bilinear_bitmap_rgb(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bilinear_bitmap_rgb.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGB src(pixelmap);

    BitMapRGB dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::bilinear(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bilinear - bitmap - rgb): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bicubic(){
    Timer<nano_t> timer;
    
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bicubic.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);

    //std::cout << "BIT DEPTH: " << img.getBitDepth() << std::endl;

    PixelMap<png_byte> dest(
        dims.height*10, 
        dims.width*10,
        dims.channels
    );
    timer.reset();
    rescale::bicubic(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bicubic): %.5f s\n", (elapsed*1e-9f));
    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bicubic_bitmap_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bicubic_bitmap_rgba.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGBA src(pixelmap);

    BitMapRGBA dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::bicubic(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bicubic - bitmap - rgba): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bicubic_bitmap_rgb(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/bicubic_bitmap_rgb.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGB src(pixelmap);

    BitMapRGB dest(dims.height*10, dims.width*10);
    timer.reset();
    rescale::bicubic(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bicubic - bitmap - rgb): %.5f s\n", (elapsed*1e-9f));

    //write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float trilinear(){
    Timer<nano_t> timer;
    
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file_l {"../../test-database/trilinear_smaller.png"};
    std::string_view to_file_L {"../../test-database/trilinear_larger.png"};
    std::string_view to_file_target {"../../test-database/trilinear_target.png"};

    png_dims dims = getDims(from_file);
    size_t height_des = dims.height*0.7; 
    size_t width_des = dims.width*0.7;

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    PixelMap<png_byte> smaller_ref;
    PixelMap<png_byte> larger_ref;
    PixelMap<png_byte> dest(height_des, width_des, dims.channels);

    timer.reset();
    rescale::rescale_by_two(smaller_ref, larger_ref, src, height_des, width_des);
    rescale::trilinear(dest, smaller_ref, larger_ref);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (trilinear): %.5f s\n", (elapsed*1e-9f));
    write_png_file(to_file_l, smaller_ref);
    write_png_file(to_file_L, larger_ref);
    write_png_file(to_file_target, dest);
    return (elapsed*1e-9f);
}

float trilinear_rgba(){
    Timer<nano_t> timer;
    
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file_l {"../../test-database/trilinear_smaller_rgba.png"};
    std::string_view to_file_L {"../../test-database/trilinear_larger_rgba.png"};
    std::string_view to_file_target {"../../test-database/trilinear_target_rgba.png"};

    png_dims dims = getDims(from_file);
    size_t height_des = dims.height*0.7; 
    size_t width_des = dims.width*0.7;

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    BitMapRGBA smaller_ref;
    BitMapRGBA larger_ref;
    BitMapRGBA dest(height_des, width_des);

    timer.reset();
    rescale::rescale_by_two(smaller_ref, larger_ref, src, height_des, width_des);
    rescale::trilinear(dest, smaller_ref, larger_ref);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (trilinear - rgba): %.5f s\n", (elapsed*1e-9f));
    write_png_file(to_file_l, smaller_ref);
    write_png_file(to_file_L, larger_ref);
    write_png_file(to_file_target, dest);
    return (elapsed*1e-9f);
}

float trilinear_rgb(){
    Timer<nano_t> timer;
    
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file_l {"../../test-database/trilinear_smaller_rgb.png"};
    std::string_view to_file_L {"../../test-database/trilinear_larger_rgb.png"};
    std::string_view to_file_target {"../../test-database/trilinear_target_rgb.png"};

    png_dims dims = getDims(from_file);
    size_t height_des = dims.height*0.7;
    size_t width_des = dims.width*0.7;

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    BitMapRGB smaller_ref;
    BitMapRGB larger_ref;
    BitMapRGB dest(height_des, width_des);

    timer.reset();
    rescale::rescale_by_two(smaller_ref, larger_ref, src, height_des, width_des);
    rescale::trilinear(dest, smaller_ref, larger_ref);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (trilinear - rgb): %.5f s\n", (elapsed*1e-9f));
    write_png_file(to_file_l, smaller_ref);
    write_png_file(to_file_L, larger_ref);
    write_png_file(to_file_target, dest);
    return (elapsed*1e-9f);
}

float difference_of_gaussians(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/diff_of_gaussians.png"};

    png_dims dims = getDims(from_file);
    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);

    PixelMap<png_byte> temp = src;
    PixelMap<png_byte> dest(src.getHeight(), src.getWidth(), src.getChannels());

    double sigma = 1;
    timer.reset();
    DoG(dest, temp, src, sigma, 2*sigma);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (diff. of gaussians): %.5f s\n", (elapsed*1e-9f));

    write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float difference_of_gaussians_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/diff_of_gaussians_rgba.png"};

    png_dims dims = getDims(from_file);
    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    
    BitMapRGBA src(pixelmap);
    BitMapRGBA temp(dims.height, dims.width);
    BitMapRGBA dest(dims.height, dims.width);
    
    double sigma = 0.2;
    timer.reset();
    DoG(dest, temp, src, sigma, 2*sigma);
    uint64_t elapsed = timer.elapsed();
    
    printf("Computation time (diff. of gaussians bitmap - rgba): %.5f s\n", (elapsed*1e-9f));
    
    write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float difference_of_gaussians_rgb(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/diff_of_gaussians_rgb.png"};

    png_dims dims = getDims(from_file);
    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    
    BitMapRGB src(pixelmap);
    BitMapRGB temp(dims.height, dims.width);
    BitMapRGB dest(dims.height, dims.width);
    
    double sigma = 0.2;
    timer.reset();
    DoG(dest, temp, src, sigma, 2*sigma);
    uint64_t elapsed = timer.elapsed();
    
    printf("Computation time (diff. of gaussians bitmap - rgb): %.5f s\n", (elapsed*1e-9f));
    
    write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

void plot_mipmap(){
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string_view to_file {"../../test-database/mip_map.png"};

    png_dims dims = getDims(from_file);
    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    //BitMapRGBA bitmap(pixelmap);

    //MipMap<BitMapRGBA> mipmap(bitmap);
    MipMap<PixelMap<png_byte>> mipmap(pixelmap);
    printf("mipmap octaves: %ld\n", mipmap.getOctaves());
    mipmap.savefig(to_file);
}

void setFileName(std::string& path, int idx){
    char fname_ext[6];
    sprintf(fname_ext, "%05d", idx);
    path += std::string(fname_ext) + ".png";
}
/*
void generate_goomer(){
    std::string_view from_file {"../../test-database/other_database/coomba.png"};
    std::string to_path {"../../test-database/mip_map"};
    png_dims dims = getDims(from_file);
    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);

    MipMap<PixelMap<png_byte>> mipmap(pixelmap);
    
    // image area
    
    size_t height_sup = mipmap[0].getHeight();
    size_t width_sup = mipmap[0].getWidth();
    size_t min_dim = std::min(height_sup, width_sup);

    for(size_t ratio = 1; ratio < min_dim; ++ratio){
        //printf("%f\n", (ratio/static_cast<float>(min_dim)));
        //mipmap.get_interpolated_img(ratio/static_cast<float>(min_dim));
    }
    //mipmap.get_interpolated_img(0);
    setFileName(to_path, 0.1);
    //std::cout << to_path << std::endl;
}
*/

#include<vector>
void speedtests(){
    /**
     * W: 1080
     * H: 1500
     * 
     * from: 800, 874 = 2 097 600
     * from: 700, 998 = 2 395 200
     * from: 836, 836 = 2 096 688
    */
    uniform_dist<uint16_t> r_height(10, 835);
    uniform_dist<uint16_t> r_width(10, 835);
    std::vector<float> buf;
    size_t iters = 100;
    buf.reserve(iters);
    for(size_t n = 0; n < iters; ++n){
        //float t_curr = bicubic(r_height(), r_width());
        printf("%ld\r", n);
        float t_curr = difference_of_gaussians_rgb();
        buf.push_back(t_curr);
    }
    float res = 0.0f;
    for(const auto &val : buf) res += val;
    printf("Mean time: %f\n", (res/iters));
}

int main(){
    /**
     * FFT
     * https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Book%3A_Signals_and_Systems_(Baraniuk_et_al.)/13%3A_Capstone_Signal_Processing_Topics/13.02%3A_The_Fast_Fourier_Transform_(FFT)
     * https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Book%3A_Electrical_Engineering_(Johnson)/05%3A_Digital_Signal_Processing/5.09%3A_Fast_Fourier_Transform_(FFT)
     * 
     * Cooley turkey algorithm
     * https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    */
    
    //nearestNeighbor();
    //nearestNeighbor_bitmap_rgba();
    //nearestNeighbor_bitmap_rgb();
    
    //bilinear();
    //bilinear_bitmap_rgba();
    //bilinear_bitmap_rgb();
    
    //bicubic();
    //bicubic_bitmap_rgba();
    //bicubic_bitmap_rgb();

    difference_of_gaussians();
    difference_of_gaussians_rgba();
    difference_of_gaussians_rgb();

    trilinear();
    trilinear_rgba();
    trilinear_rgb();

    /*
    size_t height = 5000;//2160;
    size_t width = 514;//960;
    //PixelMap<png_byte> test(height, width, 4);
    BitMapRGB test(height, width);

    uniform_dist<uint8_t> r(0, 0xff); // 0 to 255
    for(size_t i = 0; i < height; ++i){
        for(size_t j = 0; j < width; ++j){

            test(i, j) = ((r() << 0) & 0xff) |
                         ((r() << 8) & 0xff00) |
                         ((r() << 16) & 0xff0000)|
                         ((r() << 24) & 0xff000000);
            
            for(size_t c = 0; c < test.getChannels(); ++c){
                test(i, j, c) = r();
            }
        }
    }
    MipMap<PixelMap<png_byte>> mipmap(test);
    MipMap<BitMapRGB> mipmap(test);
    write_png_file("../../test-database/test.png", mipmap[0]);
    mipmap.savefig("...");
    */

    plot_mipmap();
    return 0;
}