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

void blurr(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000069_left.png"};
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
    timer.reset();
    fft(output, src, ker);
    uint64_t elapsed = timer.elapsed();

    printf("Computation time (blurr): %.5f s\n", (elapsed*1e-9f));

    printf("channels: %ld\n", output.getChannels());
    write_png_file(to_file, output);
}

void blurr_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000069_left.png"};
    std::string_view to_file {"../../test-database/blurred_rgba.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    PixelMap<png_byte> inter(src.getHeight(), src.getWidth(), src.getChannels());
    
    read_png_file(from_file, src);
    BitMapRGBA bitmap(src);

    double sd = 1.1;
    size_t k_w = std::ceil(sd * 6);
    Kernel1D<float> ker(k_w, src.getChannels(), sd);
    BitMapRGBA output(dims.height, dims.width);
    
    timer.reset();
    fft(output, bitmap, ker);
    uint64_t elapsed = timer.elapsed();

    printf("Computation time (blurr rgba): %.5f s\n", (elapsed*1e-9f));

    // output images are invisible if we take an RGB img as an RGBA img
    // so we convert the rgba to rgb to discard alpha channel
    if(dims.channels == 3){
        printf("Converting rgba image to rgb before saving\n");
        BitMapRGB bitmap_rgb(output);
        write_png_file(to_file, bitmap_rgb);
    }else{
        write_png_file(to_file, output);
    }
}

void blurr_rgb(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000069_left.png"};
    std::string_view to_file {"../../test-database/blurred_rgb.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    
    BitMapRGB bitmap(src);
    BitMapRGB rescaled(dims.height/4, dims.width/4);
    rescale::bicubic(rescaled, bitmap);

    double sd = std::pow(std::sqrt(2)/2, 1);
    size_t k_w = std::ceil(sd * 6);
    Kernel1D<float> ker(k_w, src.getChannels(), sd);

    BitMapRGB output(dims.height/4, dims.width/4);
    timer.reset();
    fft(output, rescaled, ker);
    uint64_t elapsed = timer.elapsed();

    printf("Computation time (blurr rgb): %.5f s\n", (elapsed*1e-9f));

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

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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
    write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float bicubic_bitmap_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000369_left"};
    std::string_view to_file {"../../test-database/bicubic_bitmap_rgba.png"};

    png_dims dims = getDims(from_file);

    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);
    BitMapRGBA src(pixelmap);

    BitMapRGBA dest(dims.height/2, dims.width/2);
    timer.reset();
    rescale::bicubic(dest, src);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (bicubic - bitmap - rgba): %.5f s\n", (elapsed*1e-9f));

    write_png_file(to_file, dest);
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

    write_png_file(to_file, dest);
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

    std::string_view from_file {"../../test-database/original/000069_left.png"};
    std::string_view to_file {"../../test-database/diff_of_gaussians.png"};

    png_dims dims = getDims(from_file);
    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);

    PixelMap<png_byte> temp = src;
    PixelMap<png_byte> dest(src.getHeight(), src.getWidth(), src.getChannels());

    double sigma = 0.2;
    timer.reset();
    DoG(dest, temp, src, sigma, 2*sigma);

    uint64_t elapsed = timer.elapsed();
    printf("Computation time (diff. of gaussians): %.5f s\n", (elapsed*1e-9f));

    write_png_file(to_file, dest);
    return (elapsed*1e-9f);
}

float difference_of_gaussians_rgba(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000069_left.png"};
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

    std::string_view from_file {"../../test-database/original/000069_left.png"};
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
    std::string_view from_file {"../../test-database/original/000069_left.png"};
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

inline std::string getFileName(std::string path, int idx){
    char fname_ext[6];
    sprintf(fname_ext, "%05d", idx);
    path += std::string(fname_ext) + ".png";

    return path;
}

void gen_noise(){
    size_t height = 5000;//2160;
    size_t width = 514;//960;
    //PixelMap<png_byte> test(height, width, 4);
    BitMapRGBA test(height, width);

    uniform_dist<uint8_t> r(0, 0xff); // 0 to 255
    for(size_t i = 0; i < height; ++i){
        for(size_t j = 0; j < width; ++j){
            // uncomment if using BitmapRGBA
            test(i, j) = ((r() << 0) & 0xff) |
                         ((r() << 8) & 0xff00) |
                         ((r() << 16) & 0xff0000)|
                         ((r() << 24) & 0xff000000);
            // uncomment if using PixelMap<png_byte>
            /*
            for(size_t c = 0; c < test.getChannels(); ++c){
                test(i, j, c) = r();
            }
            */
        }
    }
    //MipMap<PixelMap<png_byte>> mipmap(test);
    MipMap<BitMapRGBA> mipmap(test);
    write_png_file("../../test-database/test.png", mipmap[0]);
    mipmap.savefig("../../test-database/mipmap_noise.png");
}

void generate_goomer(){
    printf("----------------------------\n");
    std::string_view from_file {"../../test-database/original/000369_left.png"};
    std::string to_path {"../../test-database/mip_map/"};
    png_dims dims = getDims(from_file);
    PixelMap<png_byte> pixelmap(dims.height, dims.width, dims.channels);
    read_png_file(from_file, pixelmap);

    MipMap<PixelMap<png_byte>> mipmap(pixelmap);
    
    //size_t height_sup = mipmap[0].getHeight();
    //size_t width_sup = mipmap[0].getWidth();
    size_t channels = mipmap[0].getChannels();

    std::deque<PixelMap<png_byte>>::const_reverse_iterator it_r = mipmap.rbegin();
    std::deque<PixelMap<png_byte>>::const_reverse_iterator it_rend = mipmap.rend();
    
    PixelMap<png_byte> scale_target = *(it_rend-1); //(height_sup*0.5, width_sup*0.5, channels);

    size_t counter = 0;
    // keeps track of already visited dimensions (we continue loop if already visited)
    size_t mem_h = 0;
    size_t mem_w = 0;
    for(size_t i = 0; it_r != it_rend-1; ++it_r, ++i){
        size_t h = it_r->getHeight();
        size_t w = it_r->getWidth();
        size_t h_next = std::next(it_r)->getHeight();
        size_t w_next = std::next(it_r)->getWidth();
        size_t dh = h_next - h;
        size_t dw = w_next - w;
        size_t surface = dh * dw;
        
        for(size_t r = 0; r < surface; ++r){
            size_t new_h = h + r / static_cast<float>(dw);
            size_t new_w = w + r / static_cast<float>(dh);
            if(new_h == mem_h && new_w == mem_w) continue;
            else{
                std::string s = getFileName(to_path, counter);
                std::cout << "writing in: " << s << std::endl;
                PixelMap<png_byte> scale_curr(new_h, new_w, channels);
                
                rescale::trilinear(scale_curr, *it_r, *std::next(it_r));
                rescale::bicubic(scale_target, scale_curr);
                write_png_file(s, scale_target);
                ++counter;
            }
            mem_h = new_h;
            mem_w = new_w;
        }
    }
}

inline std::string get_scale_space_file(std::string path, int i, int j){
    char fname_ext[5];
    sprintf(fname_ext, "%02d%02d", i, j);
    path += std::string(fname_ext) + ".png";

    return path;
}

void generate_scale_space_rgb(){

    std::string_view from_file {"../../test-database/original/000069_left.png"};
    std::string to_folder {"../../test-database/scale_space/"};
    png_dims dims = getDims(from_file);
    PixelMap<png_byte> src(dims.height, dims.width, dims.channels);
    read_png_file(from_file, src);
    BitMapRGB bitmap(src);

    int div_h = 1;
    int div_w = 1;
    for(int i = 0; i < 4; ++i){

        double sd = std::sqrt(2)/2.0;
        for(int j = 0; j < 5; ++j){
            BitMapRGB rescaled(dims.height/div_h, dims.width/div_w);
            rescale::bicubic(rescaled, bitmap);

            size_t k_w = std::ceil(sd * 6);
            Kernel1D<float> ker(k_w, src.getChannels(), sd);
            BitMapRGB output(dims.height/div_h, dims.width/div_w);
            fft(output, rescaled, ker);

            sd *= std::sqrt(2);
            write_png_file(get_scale_space_file(to_folder, i, j), output);
        }

        div_h *= 2;
        div_w *= 2;
    }
}

void generate_DoG(){
    std::string from_folder {"../../test-database/scale_space/"};
    std::string to_folder {"../../test-database/scale_space/DoG/"};

    int div_h = 1;
    int div_w = 1;
    for(int i = 0; i < 3; ++i){

        double sd = std::sqrt(2)/2.0;
        for(int j = 0; j < 3; ++j){
            png_dims dims = getDims(get_scale_space_file(from_folder, i, j));
            
            PixelMap<png_byte> curr(dims.height, dims.width, dims.channels);
            PixelMap<png_byte> next(dims.height, dims.width, dims.channels);

            read_png_file(get_scale_space_file(from_folder, i, j), curr);
            read_png_file(get_scale_space_file(from_folder, i, j+1), next);

            BitMapRGB bitmap_curr(curr);
            BitMapRGB bitmap_next(next);

            BitMapRGB bitmap_curr_rescaled(dims.height/div_h, dims.width/div_w);
            BitMapRGB bitmap_next_rescaled(dims.height/div_h, dims.width/div_w);
            rescale::bicubic(bitmap_curr_rescaled, bitmap_curr);
            rescale::bicubic(bitmap_next_rescaled, bitmap_next);


            size_t k_w_curr = std::ceil(sd * 6);
            size_t k_w_next = std::ceil(sd * std::sqrt(2) * 6);
            Kernel1D<float> ker_curr(k_w_curr, dims.channels, sd);
            Kernel1D<float> ker_next(k_w_next, dims.channels, sd * std::sqrt(2));
            
            BitMapRGB blurred_curr(dims.height/div_h, dims.width/div_w);
            BitMapRGB blurred_next(dims.height/div_h, dims.width/div_w);

            fft(blurred_curr, bitmap_curr_rescaled, ker_curr);
            fft(blurred_next, bitmap_next_rescaled, ker_next);

            blurred_curr.subtract(blurred_next);

            sd *= std::sqrt(2);
            write_png_file(get_scale_space_file(to_folder, i, j), blurred_curr);
        }

        div_h *= 2;
        div_w *= 2;
    }
}

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

float transfer_dog_to_alpha(){
    Timer<nano_t> timer;

    std::string_view from_file {"../../test-database/original/000069_left.png"};
    std::string_view to_file {"../../test-database/diff_of_gaussians_rgba.png"};

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

    PixelMap<png_byte> res(dims.height, dims.width, 4);

    // x_new = (x - min)*(max_des - min_des)/(max-min)+min_des
    uint32_t max = 0;
    uint32_t min = 256;
    for(size_t i = 0; i < dims.height; ++i){
        for(size_t j = 0; j < dims.width; ++j){
            uint32_t avg = (dest.get(i, j, 0) + dest.get(i, j, 1) + dest.get(i, j, 2)) / 3;
            res(i, j, 3) = avg;
            max = (avg>max) ? avg : max;
            min = (avg<min) ? avg : min;
        }
    }
    for(size_t i = 0; i < dims.height; ++i){
        for(size_t j = 0; j < dims.width; ++j){
            uint32_t norm = (res(i, j, 3) - min) * 255 / (max - min);
            res(i, j, 3) = norm;
            res(i, j, 0) = norm;
            res(i, j, 1) = norm;
            res(i, j, 2) = norm;
        }
    }
    write_png_file(to_file, res);
    return (elapsed*1e-9f);
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
    /*
    nearestNeighbor();
    nearestNeighbor_bitmap_rgba();
    nearestNeighbor_bitmap_rgb();
    
    bilinear();
    bilinear_bitmap_rgba();
    bilinear_bitmap_rgb();
    
    bicubic();
    bicubic_bitmap_rgba();
    bicubic_bitmap_rgb();
    
    trilinear();
    trilinear_rgba();
    trilinear_rgb();
    */
    /*
    plot_mipmap();
    generate_goomer();
    */
    difference_of_gaussians();
    difference_of_gaussians_rgba();
    difference_of_gaussians_rgb();

    blurr();
    blurr_rgb();
    blurr_rgba();
    
    //generate_scale_space_rgb();
    //generate_DoG();
    return 0;
}