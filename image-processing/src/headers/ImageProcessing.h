#pragma once

#include "PngParser.h"
#include "Kernel.h"

#include <vector>

namespace rescale{

#define CLAMP(p, inf, sup) if(p < inf) { p = inf; } else if(p > sup) { p = sup; }

/**
 * Works well only for 2x downscaling
*/
void avg(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t new_channels = dest.getChannels();
    size_t ker_h = src.getHeight() - (new_height - 1) % src.getHeight();
    size_t ker_w = src.getWidth()  - (new_width  - 1) % src.getWidth();
    size_t surface = ker_h * ker_w;

    // size_t: each elem. of priv_buf needs enough mem. to accumulate values
    std::vector<size_t> priv_buf(new_channels, 0);

    for(size_t new_i = 0; new_i < new_height; ++new_i){
        for(size_t new_j = 0; new_j < new_width; ++new_j){
            
            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t ker_i = 0; ker_i < ker_h; ++ker_i){
                for(size_t ker_j = 0; ker_j < ker_w; ++ker_j){
                    for(size_t c = 0; c < new_channels; ++c){
                        priv_buf[c] += src(new_i+ker_i, new_j+ker_j, c);
                    }
                }
            }
            // compute avg and put in rescaled pixelmap
            for(size_t c = 0; c < new_channels; ++c){
                // after normalizing, each elem. of priv_buf returns in sizeof(T) range
                dest(new_i, new_j, c) = static_cast<png_byte>(priv_buf[c] / surface);
            }
        }
    }
}
/**
 * downscalses an image by a factor of 2
*/
void box(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t new_channels = dest.getChannels();
    size_t ker_h = 2;
    size_t ker_w = 2;
    size_t surface = ker_h * ker_w;

    // size_t: each elem. of priv_buf needs enough mem. to accumulate values
    std::vector<size_t> priv_buf(new_channels, 0);

    for(size_t new_i = 0; new_i < new_height; ++new_i){
        for(size_t new_j = 0; new_j < new_width; ++new_j){

            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t ker_i = 0; ker_i < ker_h; ++ker_i)
                for(size_t ker_j = 0; ker_j < ker_w; ++ker_j)
                    for(size_t c = 0; c < new_channels; ++c)
                        priv_buf[c] += src(ker_i+new_i*ker_h, ker_j+new_j*ker_w, c);
            // compute avg and put in rescaled pixelmap
            for(size_t c = 0; c < new_channels; ++c){
                // after normalizing, each elem. of priv_buf returns in sizeof(T) range
                dest(new_i, new_j, c) = static_cast<png_byte>(priv_buf[c] / surface);
            }
        }
    }
}

void nearestNeighbor(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t channels = dest.getChannels();

    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    // mult. by 2¹⁶ to avoid using float op. (faster)
    size_t h_ratio = (orig_height << 32) / new_h + 1;
    size_t w_ratio = (orig_width << 32) / new_w + 1;

    // precomp. lookup table that maps new image pos. to orig. image pos.
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i) h_orig_table[i] = ((i * h_ratio) >> 32);
    for(size_t j = 0; j < new_w; ++j) w_orig_table[j] = ((j * w_ratio) >> 32);

    for(size_t i = 0; i < new_h; ++i){
        for(size_t j = 0; j < new_w; ++j){
            size_t orig_i = h_orig_table[i];
            size_t orig_j = w_orig_table[j];
            #ifdef _OPENMP
                #pragma omp simd
            #endif
            for(size_t c = 0; c < channels; ++c){
                dest(i, j, c) = src(orig_i, orig_j, c);
            }
        }
    }
}

void nearestNeighbor(BitMapRGBA& dest, const BitMapRGBA& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();

    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    // mult. by 2¹⁶ to avoid using float op. (faster)
    size_t h_ratio = (orig_height << 32) / new_h + 1;
    size_t w_ratio = (orig_width << 32) / new_w + 1;

    // precomp. lookup table that maps new image pos. to orig. image pos.
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i) h_orig_table[i] = ((i * h_ratio) >> 32);
    for(size_t j = 0; j < new_w; ++j) w_orig_table[j] = ((j * w_ratio) >> 32);

    for(size_t i = 0; i < new_h; ++i){
        for(size_t j = 0; j < new_w; ++j){
            size_t orig_i = h_orig_table[i];
            size_t orig_j = w_orig_table[j];
            
            dest(i, j) = src(orig_i, orig_j);
        }
    }
}

void nearestNeighbor(BitMapRGB& dest, const BitMapRGB& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();

    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    // mult. by 2¹⁶ to avoid using float op. (faster)
    size_t h_ratio = (orig_height << 32) / new_h + 1;
    size_t w_ratio = (orig_width << 32) / new_w + 1;

    // precomp. lookup table that maps new image pos. to orig. image pos.
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i) h_orig_table[i] = ((i * h_ratio) >> 32);
    for(size_t j = 0; j < new_w; ++j) w_orig_table[j] = ((j * w_ratio) >> 32);

    for(size_t i = 0; i < new_h; ++i){
        for(size_t j = 0; j < new_w; ++j){
            size_t orig_i = h_orig_table[i];
            size_t orig_j = w_orig_table[j];
            
            dest(i, j) = src(orig_i, orig_j);
        }
    }
}

void bilinear(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();

    size_t channels = dest.getChannels();
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;
    
    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }

    // kernel mask for original height
    const ker_mask<true> km_h(2, 0, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(2, 0, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){
            
            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
            for(size_t c = 0; c < channels; ++c){

                png_byte p00 = src(orig_i  , orig_j  , c);
                png_byte p01 = src(orig_i  , orig_j+km_w(orig_j, 1), c);
                png_byte p10 = src(orig_i+km_h(orig_i, 1), orig_j  , c);
                png_byte p11 = src(orig_i+km_h(orig_i, 1), orig_j+km_w(orig_j, 1), c);

                float pixel_v = p00 * (1 - dw) * (1 - dh) + 
                                p01 *      dw  * (1 - dh) +
                                p10 *      dh  * (1 - dw) +
                                p11 *      dw  *      dh;

                dest(i, j, c) = pixel_v;
            }
        }
    }
}

void bilinear(BitMapRGBA& dest, const BitMapRGBA& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();

    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;
    
    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }

    // kernel mask for original height
    const ker_mask<true> km_h(2, 0, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(2, 0, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){
            
            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];

            uint32_t p00 = src(orig_i  , orig_j  );
            uint32_t p01 = src(orig_i  , orig_j+km_w(orig_j, 1));
            uint32_t p10 = src(orig_i+km_h(orig_i, 1), orig_j  );
            uint32_t p11 = src(orig_i+km_h(orig_i, 1), orig_j+km_w(orig_j, 1));

            uint32_t r = ( p00         & 0xff) * (1 - dw) * (1 - dh) + 
                         ( p01         & 0xff) *      dw  * (1 - dh) +
                         ( p10         & 0xff) *      dh  * (1 - dw) +
                         ( p11         & 0xff) *      dw  *      dh;

            uint32_t g = ((p00 >> 8)   & 0xff) * (1 - dw) * (1 - dh) + 
                         ((p01 >> 8)   & 0xff) *      dw  * (1 - dh) +
                         ((p10 >> 8)   & 0xff) *      dh  * (1 - dw) +
                         ((p11 >> 8)   & 0xff) *      dw  *      dh;

            uint32_t b = ((p00 >> 16)  & 0xff) * (1 - dw) * (1 - dh) + 
                         ((p01 >> 16)  & 0xff) *      dw  * (1 - dh) +
                         ((p10 >> 16)  & 0xff) *      dh  * (1 - dw) +
                         ((p11 >> 16)  & 0xff) *      dw  *      dh;
            
            uint32_t a = ((p00 >> 24)  & 0xff) * (1 - dw) * (1 - dh) + 
                         ((p01 >> 24)  & 0xff) *      dw  * (1 - dh) +
                         ((p10 >> 24)  & 0xff) *      dh  * (1 - dw) +
                         ((p11 >> 24)  & 0xff) *      dw  *      dh;
            
            dest(i, j) = ((a << 24) & 0xff000000) |
                         ((b << 16) &   0xff0000) |
                         ((g << 8 ) &     0xff00) |
                         ( r        &       0xff);
        }
    }
}

void bilinear(BitMapRGB& dest, const BitMapRGB& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();

    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;
    
    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }

    // kernel mask for original height
    const ker_mask<true> km_h(2, 0, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(2, 0, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){
            
            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];

            uint32_t p00 = src(orig_i  , orig_j  );
            uint32_t p01 = src(orig_i  , orig_j+km_w(orig_j, 1));
            uint32_t p10 = src(orig_i+km_h(orig_i, 1), orig_j  );
            uint32_t p11 = src(orig_i+km_h(orig_i, 1), orig_j+km_w(orig_j, 1));

            uint32_t r = ( p00         & 0xff) * (1 - dw) * (1 - dh) + 
                         ( p01         & 0xff) *      dw  * (1 - dh) +
                         ( p10         & 0xff) *      dh  * (1 - dw) +
                         ( p11         & 0xff) *      dw  *      dh;

            uint32_t g = ((p00 >> 8)   & 0xff) * (1 - dw) * (1 - dh) + 
                         ((p01 >> 8)   & 0xff) *      dw  * (1 - dh) +
                         ((p10 >> 8)   & 0xff) *      dh  * (1 - dw) +
                         ((p11 >> 8)   & 0xff) *      dw  *      dh;

            uint32_t b = ((p00 >> 16)  & 0xff) * (1 - dw) * (1 - dh) + 
                         ((p01 >> 16)  & 0xff) *      dw  * (1 - dh) +
                         ((p10 >> 16)  & 0xff) *      dh  * (1 - dw) +
                         ((p11 >> 16)  & 0xff) *      dw  *      dh;
            
            //CLAMP(pixel_v, 0.0f, 255.0f);
            dest(i, j) = ((b << 16) &   0xff0000) |
                         ((g << 8 ) &     0xff00) |
                         ( r        &       0xff);
        }
    }
}
/**
 * @param ref_IMG: larger reference img
 * @param ref_img: smaller reference image
*/
void trilinear(PixelMap<png_byte>& dest, 
        const PixelMap<png_byte>& ref_img, const PixelMap<png_byte>& ref_IMG){
    
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t channels = dest.getChannels();
    
    size_t HEIGHT = ref_IMG.getHeight();
    size_t WIDTH = ref_IMG.getWidth();

    size_t height = ref_img.getHeight();
    size_t width = ref_img.getWidth();
    // larger image ratio
    float H_ratio = static_cast<float>(HEIGHT-1) / new_h;
    float W_ratio = static_cast<float>(WIDTH-1) / new_w;
    // smaller image ratio
    float h_ratio = static_cast<float>(height-1) / new_h;
    float w_ratio = static_cast<float>(width-1) / new_w;
    // distance btw larger raference and interpolated image
    float dist = (WIDTH - new_w) / static_cast<float>(WIDTH - width);

    // kernel mask for original height (large img)
    const ker_mask<true> km_H(2, 0, HEIGHT);
    // kernel mask for original width (large img)
    const ker_mask<true> km_W(2, 0, WIDTH);
    // kernel mask for original height (small img)
    const ker_mask<true> km_h(2, 0, height);
    // kernel mask for original width (small img)
    const ker_mask<true> km_w(2, 0, width);
    for(size_t i = 0; i < new_h; ++i){
        // original coordinates of larger image
        size_t I_orig = H_ratio * i;
        float dH = H_ratio * i - I_orig;
        size_t i_orig = h_ratio * i;
        float dh = h_ratio * i - i_orig;
        for(size_t j = 0; j < new_w; ++j){
            // original coordinates of smaller image
            size_t J_orig = W_ratio * j;
            float dW = W_ratio * j - J_orig;
            size_t j_orig = w_ratio * j;
            float dw = w_ratio * j - j_orig;
            for(size_t c = 0; c < channels; ++c){
                // ref pixels from larger image
                png_byte P00 = ref_IMG(I_orig  , J_orig  , c);
                png_byte P01 = ref_IMG(I_orig  , J_orig+km_W(J_orig, 1), c);
                png_byte P10 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig  , c);
                png_byte P11 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig+km_W(J_orig, 1), c);
                // ref pixels from smaller image
                png_byte p00 = ref_img(i_orig  , j_orig  , c);
                png_byte p01 = ref_img(i_orig  , j_orig+km_w(j_orig, 1), c);
                png_byte p10 = ref_img(i_orig+km_h(i_orig, 1), j_orig  , c);
                png_byte p11 = ref_img(i_orig+km_h(i_orig, 1), j_orig+km_w(j_orig, 1), c);

                float pixel_v = P00 * (1 - dW) * (1 - dH) * (1 - dist) + 
                                P01 *      dW  * (1 - dH) * (1 - dist) +
                                P10 *      dH  * (1 - dW) * (1 - dist) +
                                P11 *      dW  *      dH  * (1 - dist) +
                                p00 * (1 - dw) * (1 - dh) *      dist  +
                                p01 *      dw  * (1 - dh) *      dist  +
                                p10 *      dh  * (1 - dw) *      dist  +
                                p11 *      dw  *      dh  *      dist;
                //CLAMP(pixel_v, 0.0f, 255.0f);
                dest(i, j, c) = pixel_v;
            }
        }
    }
}
/**
 * @param ref_IMG: larger reference img
 * @param ref_img: smaller reference image
*/
void trilinear(BitMapRGBA& dest, const BitMapRGBA& ref_img, const BitMapRGBA& ref_IMG){
    
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    
    size_t HEIGHT = ref_IMG.getHeight();
    size_t WIDTH = ref_IMG.getWidth();

    size_t height = ref_img.getHeight();
    size_t width = ref_img.getWidth();
    // larger image ratio
    float H_ratio = static_cast<float>(HEIGHT-1) / new_h;
    float W_ratio = static_cast<float>(WIDTH-1) / new_w;
    // smaller image ratio
    float h_ratio = static_cast<float>(height-1) / new_h;
    float w_ratio = static_cast<float>(width-1) / new_w;
    // distance btw larger raference and interpolated image
    float dist = (WIDTH - new_w) / static_cast<float>(WIDTH - width);

    // kernel mask for original height (large img)
    const ker_mask<true> km_H(2, 0, HEIGHT);
    // kernel mask for original width (large img)
    const ker_mask<true> km_W(2, 0, WIDTH);
    // kernel mask for original height (small img)
    const ker_mask<true> km_h(2, 0, height);
    // kernel mask for original width (small img)
    const ker_mask<true> km_w(2, 0, width);
    for(size_t i = 0; i < new_h; ++i){
        // larger img i
        size_t I_orig = H_ratio * i;
        float dH = H_ratio * i - I_orig;
        // smaller img i
        size_t i_orig = h_ratio * i;
        float dh = h_ratio * i - i_orig;
        for(size_t j = 0; j < new_w; ++j){
            // larger image j   
            size_t J_orig = W_ratio * j;
            float dW = W_ratio * j - J_orig;
            // smaller image j
            size_t j_orig = w_ratio * j;
            float dw = w_ratio * j - j_orig;

            // ref pixels from larger image
            uint32_t P00 = ref_IMG(I_orig  , J_orig);
            uint32_t P01 = ref_IMG(I_orig  , J_orig+km_W(J_orig, 1));
            uint32_t P10 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig);
            uint32_t P11 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig+km_W(J_orig, 1));
            // ref pixels from smaller image
            uint32_t p00 = ref_img(i_orig  , j_orig);
            uint32_t p01 = ref_img(i_orig  , j_orig+km_w(j_orig, 1));
            uint32_t p10 = ref_img(i_orig+km_h(i_orig, 1), j_orig);
            uint32_t p11 = ref_img(i_orig+km_h(i_orig, 1), j_orig+km_w(j_orig, 1));

            uint32_t r = ( P00         & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ( P01         & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ( P10         & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ( P11         & 0xff) *      dW  *      dH  * (1 - dist) +

                         ( p00         & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ( p01         & 0xff) *      dw  * (1 - dh) *      dist  +
                         ( p10         & 0xff) *      dh  * (1 - dw) *      dist  +
                         ( p11         & 0xff) *      dw  *      dh  *      dist;

            uint32_t g = ((P00 >> 8)   & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ((P01 >> 8)   & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ((P10 >> 8)   & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ((P11 >> 8)   & 0xff) *      dW  *      dH  * (1 - dist) +

                         ((p00 >> 8)   & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ((p01 >> 8)   & 0xff) *      dw  * (1 - dh) *      dist  +
                         ((p10 >> 8)   & 0xff) *      dh  * (1 - dw) *      dist  +
                         ((p11 >> 8)   & 0xff) *      dw  *      dh  *      dist;

            uint32_t b = ((P00 >> 16)  & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ((P01 >> 16)  & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ((P10 >> 16)  & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ((P11 >> 16)  & 0xff) *      dW  *      dH  * (1 - dist) +

                         ((p00 >> 16)  & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ((p01 >> 16)  & 0xff) *      dw  * (1 - dh) *      dist  +
                         ((p10 >> 16)  & 0xff) *      dh  * (1 - dw) *      dist  +
                         ((p11 >> 16)  & 0xff) *      dw  *      dh  *      dist;
            
            uint32_t a = ((P00 >> 24)  & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ((P01 >> 24)  & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ((P10 >> 24)  & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ((P11 >> 24)  & 0xff) *      dW  *      dH  * (1 - dist) +

                         ((p00 >> 24)  & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ((p01 >> 24)  & 0xff) *      dw  * (1 - dh) *      dist  +
                         ((p10 >> 24)  & 0xff) *      dh  * (1 - dw) *      dist  +
                         ((p11 >> 24)  & 0xff) *      dw  *      dh  *      dist;
            
            //CLAMP(pixel_v, 0.0f, 255.0f);
            dest(i, j) = ((a << 24) & 0xff000000) |
                         ((b << 16) &   0xff0000) |
                         ((g << 8 ) &     0xff00) |
                         ( r        &       0xff);
        }
    }
}
/**
 * @param ref_IMG: larger reference img
 * @param ref_img: smaller reference image
*/
void trilinear(BitMapRGB& dest, const BitMapRGB& ref_img, const BitMapRGB& ref_IMG){
    
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    
    size_t HEIGHT = ref_IMG.getHeight();
    size_t WIDTH = ref_IMG.getWidth();

    size_t height = ref_img.getHeight();
    size_t width = ref_img.getWidth();
    // larger image ratio
    float H_ratio = static_cast<float>(HEIGHT-1) / new_h;
    float W_ratio = static_cast<float>(WIDTH-1) / new_w;
    // smaller image ratio
    float h_ratio = static_cast<float>(height-1) / new_h;
    float w_ratio = static_cast<float>(width-1) / new_w;
    // distance btw larger raference and interpolated image
    float dist = (WIDTH - new_w) / static_cast<float>(WIDTH - width);

    // kernel mask for original height (large img)
    const ker_mask<true> km_H(2, 0, HEIGHT);
    // kernel mask for original width (large img)
    const ker_mask<true> km_W(2, 0, WIDTH);
    // kernel mask for original height (small img)
    const ker_mask<true> km_h(2, 0, height);
    // kernel mask for original width (small img)
    const ker_mask<true> km_w(2, 0, width);

    for(size_t i = 0; i < new_h; ++i){
        // larger img i
        size_t I_orig = H_ratio * i;
        float dH = H_ratio * i - I_orig;
        // smaller img i
        size_t i_orig = h_ratio * i;
        float dh = h_ratio * i - i_orig;
        
        for(size_t j = 0; j < new_w; ++j){
            // larger image j   
            size_t J_orig = W_ratio * j;
            float dW = W_ratio * j - J_orig;
            // smaller image j
            size_t j_orig = w_ratio * j;
            float dw = w_ratio * j - j_orig;

            // ref pixels from larger image
            uint32_t P00 = ref_IMG(I_orig  , J_orig);
            uint32_t P01 = ref_IMG(I_orig  , J_orig+km_W(J_orig, 1));
            uint32_t P10 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig);
            uint32_t P11 = ref_IMG(I_orig+km_H(I_orig, 1), J_orig+km_W(J_orig, 1));
            // ref pixels from smaller image
            uint32_t p00 = ref_img(i_orig  , j_orig);
            uint32_t p01 = ref_img(i_orig  , j_orig+km_w(j_orig, 1));
            uint32_t p10 = ref_img(i_orig+km_h(i_orig, 1), j_orig);
            uint32_t p11 = ref_img(i_orig+km_h(i_orig, 1), j_orig+km_w(j_orig, 1));

            uint32_t r = ( P00         & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ( P01         & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ( P10         & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ( P11         & 0xff) *      dW  *      dH  * (1 - dist) +

                         ( p00         & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ( p01         & 0xff) *      dw  * (1 - dh) *      dist  +
                         ( p10         & 0xff) *      dh  * (1 - dw) *      dist  +
                         ( p11         & 0xff) *      dw  *      dh  *      dist;

            uint32_t g = ((P00 >> 8)   & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ((P01 >> 8)   & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ((P10 >> 8)   & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ((P11 >> 8)   & 0xff) *      dW  *      dH  * (1 - dist) +

                         ((p00 >> 8)   & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ((p01 >> 8)   & 0xff) *      dw  * (1 - dh) *      dist  +
                         ((p10 >> 8)   & 0xff) *      dh  * (1 - dw) *      dist  +
                         ((p11 >> 8)   & 0xff) *      dw  *      dh  *      dist;

            uint32_t b = ((P00 >> 16)  & 0xff) * (1 - dW) * (1 - dH) * (1 - dist) + 
                         ((P01 >> 16)  & 0xff) *      dW  * (1 - dH) * (1 - dist) + 
                         ((P10 >> 16)  & 0xff) *      dH  * (1 - dW) * (1 - dist) + 
                         ((P11 >> 16)  & 0xff) *      dW  *      dH  * (1 - dist) +

                         ((p00 >> 16)  & 0xff) * (1 - dw) * (1 - dh) *      dist  + 
                         ((p01 >> 16)  & 0xff) *      dw  * (1 - dh) *      dist  +
                         ((p10 >> 16)  & 0xff) *      dh  * (1 - dw) *      dist  +
                         ((p11 >> 16)  & 0xff) *      dw  *      dh  *      dist;
            
            //CLAMP(pixel_v, 0.0f, 255.0f);
            dest(i, j) = ((b << 16) &   0xff0000) |
                         ((g << 8 ) &     0xff00) |
                         ( r        &       0xff);
        }
    }
}

inline float cubicHermite(const float p0, const float p1, const float p2, const float p3, const float val){
    return p1 + 0.5f*val*(p2 - p0 + val*(2.0f*p0 - 5.0f*p1 + 4.0f*p2 - p3 + val*(3.0f*(p1 - p2) + p3 - p0)));
}

void bicubic(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;

    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }
    // kernel mask for original height
    const ker_mask<true> km_h(4, 1, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(4, 1, orig_width);
    for(size_t i = 0; i < new_h; ++i){
        
        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        //size_t orig_i = i * h_ratio;
        //float dh = h_ratio * i - orig_i;
        for(size_t j = 0; j < new_w; ++j){

            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
            //size_t orig_j = j * w_ratio;
            //float dw = w_ratio * j - orig_j;

            const png_byte* p00 = src.pixBegin(orig_i-km_h(orig_i, 0), orig_j-km_w(orig_i, 0));
            const png_byte* p01 = src.pixBegin(orig_i-km_h(orig_i, 0), orig_j  );
            const png_byte* p02 = src.pixBegin(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_i, 2));
            const png_byte* p03 = src.pixBegin(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_i, 3));

            const png_byte* p10 = src.pixBegin(orig_i  , orig_j-km_w(orig_i, 0));
            const png_byte* p11 = src.pixBegin(orig_i  , orig_j  );
            const png_byte* p12 = src.pixBegin(orig_i  , orig_j+km_w(orig_i, 2));
            const png_byte* p13 = src.pixBegin(orig_i  , orig_j+km_w(orig_i, 3));
                
            const png_byte* p20 = src.pixBegin(orig_i+km_h(orig_i, 2), orig_j-km_w(orig_i, 0));
            const png_byte* p21 = src.pixBegin(orig_i+km_h(orig_i, 2), orig_j  );
            const png_byte* p22 = src.pixBegin(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_i, 2));
            const png_byte* p23 = src.pixBegin(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_i, 3));
                
            const png_byte* p30 = src.pixBegin(orig_i+km_h(orig_i, 3), orig_j-km_w(orig_i, 0));
            const png_byte* p31 = src.pixBegin(orig_i+km_h(orig_i, 3), orig_j  );
            const png_byte* p32 = src.pixBegin(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_i, 2));
            const png_byte* p33 = src.pixBegin(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_i, 3));

            png_byte* pixel_dest = dest.pixBegin(i, j);

            for(; pixel_dest != dest.pixEnd(i, j);  ++p00, ++p01, ++p02, ++p03, 
                                                    ++p10, ++p11, ++p12, ++p13, 
                                                    ++p20, ++p21, ++p22, ++p23,
                                                    ++p30, ++p31, ++p32, ++p33,
                                                    ++pixel_dest){
                const float p0 = cubicHermite(*p00, *p01, *p02, *p03, dw);
                const float p1 = cubicHermite(*p10, *p11, *p12, *p13, dw);
                const float p2 = cubicHermite(*p20, *p21, *p22, *p23, dw);
                const float p3 = cubicHermite(*p30, *p31, *p32, *p33, dw);
                
                float pixel_v_c = cubicHermite(p0, p1, p2, p3, dh);
                    
                CLAMP(pixel_v_c, 0.0f, 255.0f);
                    
                *pixel_dest = pixel_v_c;
            }
        }
    }
}

void bicubic_loop(BitMapRGBA& dest, const BitMapRGBA& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;

    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);

    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }

    // kernel mask for original height
    const ker_mask<true> km_h(4, 1, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(4, 1, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){

            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
                
            const uint32_t p00 = src(orig_i-km_h(orig_i, 0), orig_j-km_w(orig_i, 0));
            const uint32_t p01 = src(orig_i-km_h(orig_i, 0), orig_j  );
            const uint32_t p02 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_i, 2));
            const uint32_t p03 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_i, 3));

            const uint32_t p10 = src(orig_i  , orig_j-km_w(orig_i, 0));
            const uint32_t p11 = src(orig_i  , orig_j  );
            const uint32_t p12 = src(orig_i  , orig_j+km_w(orig_i, 2));
            const uint32_t p13 = src(orig_i  , orig_j+km_w(orig_i, 3));
                
            const uint32_t p20 = src(orig_i+km_h(orig_i, 2), orig_j-km_w(orig_i, 0));
            const uint32_t p21 = src(orig_i+km_h(orig_i, 2), orig_j  );
            const uint32_t p22 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_i, 2));
            const uint32_t p23 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_i, 3));
                
            const uint32_t p30 = src(orig_i+km_h(orig_i, 3), orig_j-km_w(orig_i, 0));
            const uint32_t p31 = src(orig_i+km_h(orig_i, 3), orig_j  );
            const uint32_t p32 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_i, 2));
            const uint32_t p33 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_i, 3));

            dest(i, j) = 0; // making sure pixel is zero to start accumulating
            for(size_t c = 0; c < 4; ++c){
                const float p0 = cubicHermite((p00 >> 8*c & 0xff), (p01 >> 8*c & 0xff), (p02 >> 8*c & 0xff), (p03 >> 8*c & 0xff), dw);
                const float p1 = cubicHermite((p10 >> 8*c & 0xff), (p11 >> 8*c & 0xff), (p12 >> 8*c & 0xff), (p13 >> 8*c & 0xff), dw);
                const float p2 = cubicHermite((p20 >> 8*c & 0xff), (p21 >> 8*c & 0xff), (p22 >> 8*c & 0xff), (p23 >> 8*c & 0xff), dw);
                const float p3 = cubicHermite((p30 >> 8*c & 0xff), (p31 >> 8*c & 0xff), (p32 >> 8*c & 0xff), (p33 >> 8*c & 0xff), dw);

                float p = cubicHermite(p0, p1, p2, p3, dh);
                CLAMP(p, 0.0f, 255.0f);

                dest(i, j) |= ((static_cast<uint32_t>(p) << rgba_utils.shift[c]) & rgba_utils.rgba[c]);
            }
        }
    }
}

void bicubic(BitMapRGBA& dest, const BitMapRGBA& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;

    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);

    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }
    // kernel mask for original height
    const ker_mask<true> km_h(4, 1, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(4, 1, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){
            
            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
                
            const uint32_t p00 = src(orig_i-km_h(orig_i, 0), orig_j-km_w(orig_j, 0));
            const uint32_t p01 = src(orig_i-km_h(orig_i, 0), orig_j  );
            const uint32_t p02 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 2));
            const uint32_t p03 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 3));

            const uint32_t p10 = src(orig_i  , orig_j-km_w(orig_j, 0));
            const uint32_t p11 = src(orig_i  , orig_j  );
            const uint32_t p12 = src(orig_i  , orig_j+km_w(orig_j, 2));
            const uint32_t p13 = src(orig_i  , orig_j+km_w(orig_j, 3));
                
            const uint32_t p20 = src(orig_i+km_h(orig_i, 2), orig_j-km_w(orig_j, 0));
            const uint32_t p21 = src(orig_i+km_h(orig_i, 2), orig_j  );
            const uint32_t p22 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 2));
            const uint32_t p23 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 3));
                
            const uint32_t p30 = src(orig_i+km_h(orig_i, 3), orig_j-km_w(orig_j, 0));
            const uint32_t p31 = src(orig_i+km_h(orig_i, 3), orig_j  );
            const uint32_t p32 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 2));
            const uint32_t p33 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 3));
            
            const float r0 = cubicHermite((p00 >> 0 & 0xff), (p01 >> 0 & 0xff), (p02 >> 0 & 0xff), (p03 >> 0 & 0xff), dw);
            const float r1 = cubicHermite((p10 >> 0 & 0xff), (p11 >> 0 & 0xff), (p12 >> 0 & 0xff), (p13 >> 0 & 0xff), dw);
            const float r2 = cubicHermite((p20 >> 0 & 0xff), (p21 >> 0 & 0xff), (p22 >> 0 & 0xff), (p23 >> 0 & 0xff), dw);
            const float r3 = cubicHermite((p30 >> 0 & 0xff), (p31 >> 0 & 0xff), (p32 >> 0 & 0xff), (p33 >> 0 & 0xff), dw);

            const float g0 = cubicHermite((p00 >> 8 & 0xff), (p01 >> 8 & 0xff), (p02 >> 8 & 0xff), (p03 >> 8 & 0xff), dw);
            const float g1 = cubicHermite((p10 >> 8 & 0xff), (p11 >> 8 & 0xff), (p12 >> 8 & 0xff), (p13 >> 8 & 0xff), dw);
            const float g2 = cubicHermite((p20 >> 8 & 0xff), (p21 >> 8 & 0xff), (p22 >> 8 & 0xff), (p23 >> 8 & 0xff), dw);
            const float g3 = cubicHermite((p30 >> 8 & 0xff), (p31 >> 8 & 0xff), (p32 >> 8 & 0xff), (p33 >> 8 & 0xff), dw);

            const float b0 = cubicHermite((p00 >> 16 & 0xff), (p01 >> 16 & 0xff), (p02 >> 16 & 0xff), (p03 >> 16 & 0xff), dw);
            const float b1 = cubicHermite((p10 >> 16 & 0xff), (p11 >> 16 & 0xff), (p12 >> 16 & 0xff), (p13 >> 16 & 0xff), dw);
            const float b2 = cubicHermite((p20 >> 16 & 0xff), (p21 >> 16 & 0xff), (p22 >> 16 & 0xff), (p23 >> 16 & 0xff), dw);
            const float b3 = cubicHermite((p30 >> 16 & 0xff), (p31 >> 16 & 0xff), (p32 >> 16 & 0xff), (p33 >> 16 & 0xff), dw);

            const float a0 = cubicHermite((p00 >> 24 & 0xff), (p01 >> 24 & 0xff), (p02 >> 24 & 0xff), (p03 >> 24 & 0xff), dw);
            const float a1 = cubicHermite((p10 >> 24 & 0xff), (p11 >> 24 & 0xff), (p12 >> 24 & 0xff), (p13 >> 24 & 0xff), dw);
            const float a2 = cubicHermite((p20 >> 24 & 0xff), (p21 >> 24 & 0xff), (p22 >> 24 & 0xff), (p23 >> 24 & 0xff), dw);
            const float a3 = cubicHermite((p30 >> 24 & 0xff), (p31 >> 24 & 0xff), (p32 >> 24 & 0xff), (p33 >> 24 & 0xff), dw);

            float r = cubicHermite(r0, r1, r2, r3, dh);
            float g = cubicHermite(g0, g1, g2, g3, dh);
            float b = cubicHermite(b0, b1, b2, b3, dh);
            float a = cubicHermite(a0, a1, a2, a3, dh);

            CLAMP(r, 0.0f, 255.0f);
            CLAMP(g, 0.0f, 255.0f);
            CLAMP(b, 0.0f, 255.0f);
            CLAMP(a, 0.0f, 255.0f);
            
            dest(i, j) = ((static_cast<uint32_t>(a) << 24) & 0xff000000) |
                         ((static_cast<uint32_t>(b) << 16) &   0xff0000) |
                         ((static_cast<uint32_t>(g) << 8 ) &     0xff00) |
                         ( static_cast<uint32_t>(r)        &       0xff);
        }
    }
}

void bicubic(BitMapRGB& dest, const BitMapRGB& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;

    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
    }

    // kernel mask for original height
    const ker_mask<true> km_h(4, 1, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(4, 1, orig_width);
    for(size_t i = 0; i < new_h; ++i){

        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){
            
            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
                
            const uint32_t p00 = src(orig_i-km_h(orig_i, 0), orig_j-km_w(orig_j, 0));
            const uint32_t p01 = src(orig_i-km_h(orig_i, 0), orig_j  );
            const uint32_t p02 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 2));
            const uint32_t p03 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 3));

            const uint32_t p10 = src(orig_i  , orig_j-km_w(orig_j, 0));
            const uint32_t p11 = src(orig_i  , orig_j  );
            const uint32_t p12 = src(orig_i  , orig_j+km_w(orig_j, 2));
            const uint32_t p13 = src(orig_i  , orig_j+km_w(orig_j, 3));
                
            const uint32_t p20 = src(orig_i+km_h(orig_i, 2), orig_j-km_w(orig_j, 0));
            const uint32_t p21 = src(orig_i+km_h(orig_i, 2), orig_j  );
            const uint32_t p22 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 2));
            const uint32_t p23 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 3));
                
            const uint32_t p30 = src(orig_i+km_h(orig_i, 3), orig_j-km_w(orig_j, 0));
            const uint32_t p31 = src(orig_i+km_h(orig_i, 3), orig_j  );
            const uint32_t p32 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 2));
            const uint32_t p33 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 3));
            
            const float r0 = cubicHermite((p00 >> 0 & 0xff), (p01 >> 0 & 0xff), (p02 >> 0 & 0xff), (p03 >> 0 & 0xff), dw);
            const float r1 = cubicHermite((p10 >> 0 & 0xff), (p11 >> 0 & 0xff), (p12 >> 0 & 0xff), (p13 >> 0 & 0xff), dw);
            const float r2 = cubicHermite((p20 >> 0 & 0xff), (p21 >> 0 & 0xff), (p22 >> 0 & 0xff), (p23 >> 0 & 0xff), dw);
            const float r3 = cubicHermite((p30 >> 0 & 0xff), (p31 >> 0 & 0xff), (p32 >> 0 & 0xff), (p33 >> 0 & 0xff), dw);

            const float g0 = cubicHermite((p00 >> 8 & 0xff), (p01 >> 8 & 0xff), (p02 >> 8 & 0xff), (p03 >> 8 & 0xff), dw);
            const float g1 = cubicHermite((p10 >> 8 & 0xff), (p11 >> 8 & 0xff), (p12 >> 8 & 0xff), (p13 >> 8 & 0xff), dw);
            const float g2 = cubicHermite((p20 >> 8 & 0xff), (p21 >> 8 & 0xff), (p22 >> 8 & 0xff), (p23 >> 8 & 0xff), dw);
            const float g3 = cubicHermite((p30 >> 8 & 0xff), (p31 >> 8 & 0xff), (p32 >> 8 & 0xff), (p33 >> 8 & 0xff), dw);

            const float b0 = cubicHermite((p00 >> 16 & 0xff), (p01 >> 16 & 0xff), (p02 >> 16 & 0xff), (p03 >> 16 & 0xff), dw);
            const float b1 = cubicHermite((p10 >> 16 & 0xff), (p11 >> 16 & 0xff), (p12 >> 16 & 0xff), (p13 >> 16 & 0xff), dw);
            const float b2 = cubicHermite((p20 >> 16 & 0xff), (p21 >> 16 & 0xff), (p22 >> 16 & 0xff), (p23 >> 16 & 0xff), dw);
            const float b3 = cubicHermite((p30 >> 16 & 0xff), (p31 >> 16 & 0xff), (p32 >> 16 & 0xff), (p33 >> 16 & 0xff), dw);

            float r = cubicHermite(r0, r1, r2, r3, dh);
            float g = cubicHermite(g0, g1, g2, g3, dh);
            float b = cubicHermite(b0, b1, b2, b3, dh);

            CLAMP(r, 0.0f, 255.0f);
            CLAMP(g, 0.0f, 255.0f);
            CLAMP(b, 0.0f, 255.0f);
            
            dest(i, j) = ((static_cast<uint32_t>(b) << 16) &   0xff0000) |
                         ((static_cast<uint32_t>(g) << 8 ) &     0xff00) |
                         ( static_cast<uint32_t>(r)        &       0xff);
        }
    }
}

/**
 * This version is much less efficient than bicubic() for a single call
*/
void bicubic_cache(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src){
    size_t new_h = dest.getHeight();
    size_t new_w = dest.getWidth();
    size_t channels = dest.getChannels();
    
    size_t orig_height = src.getHeight();
    size_t orig_width = src.getWidth();
    
    float h_ratio = static_cast<float>(orig_height-1) / new_h;
    float w_ratio = static_cast<float>(orig_width-1) / new_w;
    
    png_byte p00, p01, p02, p03;
    png_byte p10, p11, p12, p13;
    png_byte p20, p21, p22, p23;
    png_byte p30, p31, p32, p33;
    
    float a00, a01, a02, a03;
    float a10, a11, a12, a13;
    float a20, a21, a22, a23;
    float a30, a31, a32, a33;

    // precompute values (a little bit faster)
    std::vector<size_t> h_orig_table(new_h);
    std::vector<size_t> w_orig_table(new_w);
    std::vector<float> dh_orig_table(new_h);
    std::vector<float> dw_orig_table(new_w);
    std::vector<float> dh2_orig_table(new_h);
    std::vector<float> dw2_orig_table(new_w);
    std::vector<float> dh3_orig_table(new_h);
    std::vector<float> dw3_orig_table(new_w);
    for(size_t i = 0; i < new_h; ++i){
        h_orig_table[i] = i * h_ratio;
        dh_orig_table[i] = h_ratio * i - h_orig_table[i];
        dh2_orig_table[i] = dh_orig_table[i] * dh_orig_table[i];
        dh3_orig_table[i] = dh2_orig_table[i] * dh_orig_table[i];
    }
    for(size_t j = 0; j < new_w; ++j){
        w_orig_table[j] = j * w_ratio;
        dw_orig_table[j] = w_ratio * j - w_orig_table[j];
        dw2_orig_table[j] = dw_orig_table[j] * dw_orig_table[j];
        dw3_orig_table[j] = dw2_orig_table[j] * dw_orig_table[j];
    }
    
    float pixel_v_c;

    // kernel mask for original height
    const ker_mask<true> km_h(4, 1, orig_height);
    // kernel mask for original width
    const ker_mask<true> km_w(4, 1, orig_width);
    for(size_t i = 0; i < new_h; ++i){
        
        size_t orig_i = h_orig_table[i];
        float dh = dh_orig_table[i];
        float dh2 = dh2_orig_table[i];
        float dh3 = dh3_orig_table[i];
        for(size_t j = 0; j < new_w; ++j){

            size_t orig_j = w_orig_table[j];
            float dw = dw_orig_table[j];
            float dw2 = dw2_orig_table[j];
            float dw3 = dw3_orig_table[j];
            for(size_t c = 0; c < channels; ++c){

                p00 = src(orig_i-km_h(orig_i, 0), orig_j-km_w(orig_j, 0), c);
                p01 = src(orig_i-km_h(orig_i, 0), orig_j  , c);
                p02 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 2), c);
                p03 = src(orig_i-km_h(orig_i, 0), orig_j+km_w(orig_j, 3), c);

                p10 = src(orig_i  , orig_j-km_w(orig_j, 0), c);
                p11 = src(orig_i  , orig_j  , c);
                p12 = src(orig_i  , orig_j+km_w(orig_j, 2), c);
                p13 = src(orig_i  , orig_j+km_w(orig_j, 3), c);
                
                p20 = src(orig_i+km_h(orig_i, 2), orig_j-km_w(orig_j, 0), c);
                p21 = src(orig_i+km_h(orig_i, 2), orig_j  , c);
                p22 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 2), c);
                p23 = src(orig_i+km_h(orig_i, 2), orig_j+km_w(orig_j, 3), c);
                
                p30 = src(orig_i+km_h(orig_i, 3), orig_j-km_w(orig_j, 0), c);
                p31 = src(orig_i+km_h(orig_i, 3), orig_j  , c);
                p32 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 2), c);
                p33 = src(orig_i+km_h(orig_i, 3), orig_j+km_w(orig_j, 3), c);

                a00 = p11;
                a01 = -.5*p10 + .5*p12;
                a02 = p10 - 2.5 * p11 + 2 * p12 - .5 * p13;
                a03 = -.5*p10 + 1.5*p11 - 1.5*p12 + .5*p13;
                a10 = -.5*p01 + .5*p21;
                a11 = .25*p00 - .25*p02 - .25*p20 + .25*p22;
                a12 = -.5*p00 + 1.25*p01 - p02 + .25*p03 + .5*p20 - 1.25*p21 + p22 - .25*p23;
                a13 = .25*p00 - .75*p01 + .75*p02 - .25*p03 - .25*p20 + .75*p21 - .75*p22 + .25*p23;
                a20 = p01 - 2.5*p11 + 2*p21 - .5*p31;
                a21 = -.5*p00 + .5*p02 + 1.25*p10 - 1.25*p12 - p20 + p22 + .25*p30 - .25*p32;
                a22 = p00 - 2.5*p01 + 2*p02 - .5*p03 - 2.5*p10 + 6.25*p11 - 5*p12 + 1.25*p13 + 2*p20 - 5*p21 + 4*p22 - p23 - .5*p30 + 1.25*p31 - p32 + .25*p33;
                a23 = -.5*p00 + 1.5*p01 - 1.5*p02 + .5*p03 + 1.25*p10 - 3.75*p11 + 3.75*p12 - 1.25*p13 - p20 + 3*p21 - 3*p22 + p23 + .25*p30 - .75*p31 + .75*p32 - .25*p33;
                a30 = -.5*p01 + 1.5*p11 - 1.5*p21 + .5*p31;
                a31 = .25*p00 - .25*p02 - .75*p10 + .75*p12 + .75*p20 - .75*p22 - .25*p30 + .25*p32;
                a32 = -.5*p00 + 1.25*p01 - p02 + .25*p03 + 1.5*p10 - 3.75*p11 + 3*p12 - .75*p13 - 1.5*p20 + 3.75*p21 - 3*p22 + .75*p23 + .5*p30 - 1.25*p31 + p32 - .25*p33;
                a33 = .25*p00 - .75*p01 + .75*p02 - .25*p03 - .75*p10 + 2.25*p11 - 2.25*p12 + .75*p13 + .75*p20 - 2.25*p21 + 2.25*p22 - .75*p23 - .25*p30 + .75*p31 - .75*p32 + .25*p33;

                pixel_v_c = (a00 + a01*dw + a02*dw2 + a03*dw3)       +
                            (a10 + a11*dw + a12*dw2 + a13*dw3) * dh  +
                            (a20 + a21*dw + a22*dw2 + a23*dw3) * dh2 +
                            (a30 + a31*dw + a32*dw2 + a33*dw3) * dh3; 
                CLAMP(pixel_v_c, 0.0f, 255.0f);
                dest(i, j, c) = pixel_v_c;
                }
            }
        }
}
/**
 * [down/up]scales by a factor of 2 until target image dims
 * are > smaller reference image dims and
 * are < larger reference image dims
 * @param h_des: desired height
 * @param w_des: desired width
 * @return ratio: how much it has been [down/up]scaled
*/
std::pair<float, float> 
rescale_by_two(PixelMap<png_byte>& smaller_ref, PixelMap<png_byte>& larger_ref, 
        const PixelMap<png_byte>& src, size_t h_des, size_t w_des){
    
    size_t h_curr = src.getHeight();
    size_t w_curr = src.getWidth();
    float w_fact = (w_des > w_curr) ? 2.0f : 0.5f;
    float h_fact = (h_des > h_curr) ? 2.0f : 0.5f;
    float w_inf = w_curr;
    float h_inf = h_curr;
    float w_sup = 2 * w_inf;
    float h_sup = 2 * h_inf;
    float w_ratio = w_fact;
    float h_ratio = h_fact;

    PixelMap<png_byte> res = src;
    while(!(h_inf < h_des && h_des <= h_sup)){
        h_inf *= h_fact;
        h_sup = 2 * h_inf;
        h_ratio *= h_fact;
        auto temp = PixelMap<png_byte>(res.getHeight()*h_fact, res.getWidth(), res.getChannels());
        bilinear(temp, res);
        res = temp;
    }
    while(!(w_inf < w_des && w_des <= w_sup)){
        w_inf *= w_fact;
        w_sup = 2 * w_inf;
        w_ratio *= w_fact;
        auto temp = PixelMap<png_byte>(res.getHeight(), res.getWidth()*w_fact, res.getChannels());
        bilinear(temp, res);
        res = temp;
    }
    smaller_ref = res;
    larger_ref = PixelMap<png_byte>(smaller_ref.getHeight()*2, smaller_ref.getWidth()*2, smaller_ref.getChannels());
    bilinear(larger_ref, smaller_ref);
    //printf("%ld, %ld, %ld\n", inf, w_des, sup);
    return std::make_pair(h_ratio, w_ratio);
}
/**
 * [down/up]scales by a factor of 2 until target image dims
 * are > smaller reference image dims and
 * are < larger reference image dims
 * @param h_des: desired height
 * @param w_des: desired width
 * @return ratio: how much it has been [down/up]scaled
*/
std::pair<float, float> 
rescale_by_two(BitMapRGBA& smaller_ref, BitMapRGBA& larger_ref, 
        const BitMapRGBA& src, size_t h_des, size_t w_des){
    
    size_t h_curr = src.getHeight();
    size_t w_curr = src.getWidth();
    float w_fact = (w_des > w_curr) ? 2.0f : 0.5f;
    float h_fact = (h_des > h_curr) ? 2.0f : 0.5f;
    float w_inf = w_curr;
    float h_inf = h_curr;
    float w_sup = 2 * w_inf;
    float h_sup = 2 * h_inf;
    float w_ratio = w_fact;
    float h_ratio = h_fact;

    BitMapRGBA res = src;
    while(!(h_inf < h_des && h_des <= h_sup)){
        h_inf *= h_fact;
        h_sup = 2 * h_inf;
        h_ratio *= h_fact;
        BitMapRGBA temp(res.getHeight()*h_fact, res.getWidth());
        bilinear(temp, res);
        res = temp;
    }
    while(!(w_inf < w_des && w_des <= w_sup)){
        w_inf *= w_fact;
        w_sup = 2 * w_inf;
        w_ratio *= w_fact;
        BitMapRGBA temp(res.getHeight(), res.getWidth()*w_fact);
        bilinear(temp, res);
        res = temp;
    }
    smaller_ref = res;
    larger_ref = BitMapRGBA(smaller_ref.getHeight()*2, smaller_ref.getWidth()*2);
    bilinear(larger_ref, smaller_ref);
    //printf("%ld, %ld, %ld\n", inf, w_des, sup);
    return std::make_pair(h_ratio, w_ratio);
}
    /**
 * [down/up]scales by a factor of 2 until target image dims
 * are > smaller reference image dims and
 * are < larger reference image dims
 * @param h_des: desired height
 * @param w_des: desired width
 * @return ratio: how much it has been [down/up]scaled
*/
std::pair<float, float> 
rescale_by_two(BitMapRGB& smaller_ref, BitMapRGB& larger_ref, 
        const BitMapRGB& src, size_t h_des, size_t w_des){
    
    size_t h_curr = src.getHeight();
    size_t w_curr = src.getWidth();
    float w_fact = (w_des > w_curr) ? 2.0f : 0.5f;
    float h_fact = (h_des > h_curr) ? 2.0f : 0.5f;
    float w_inf = w_curr;
    float h_inf = h_curr;
    float w_sup = 2 * w_inf;
    float h_sup = 2 * h_inf;
    float w_ratio = w_fact;
    float h_ratio = h_fact;

    BitMapRGB res = src;
    while(!(h_inf < h_des && h_des <= h_sup)){
        h_inf *= h_fact;
        h_sup = 2 * h_inf;
        h_ratio *= h_fact;
        BitMapRGB temp(res.getHeight()*h_fact, res.getWidth());
        bilinear(temp, res);
        res = temp;
    }
    while(!(w_inf < w_des && w_des <= w_sup)){
        w_inf *= w_fact;
        w_sup = 2 * w_inf;
        w_ratio *= w_fact;
        BitMapRGB temp(res.getHeight(), res.getWidth()*w_fact);
        bilinear(temp, res);
        res = temp;
    }
    smaller_ref = res;
    larger_ref = BitMapRGB(smaller_ref.getHeight()*2, smaller_ref.getWidth()*2);
    bilinear(larger_ref, smaller_ref);
    //printf("%ld, %ld, %ld\n", inf, w_des, sup);
    return std::make_pair(h_ratio, w_ratio);
}
}
/**
 * The kernels for blurring are separable
 * 1 kernel slides [horizont/vertic]aly
 * 
 * OUT: output tensor
 * IN: intput tensor
 * K: kernel
 * k: kernel length
 * w: input/output width (same for blurr)
 * 
 * case: kernel length ODD (ker. starts with 1 pix. out and ends with 1 pix. out)
 *      OUT[0     ->   k/2) = IN[0     ->   k/2) * K[  k/2 -> k  )
 *      OUT[k/2   -> w-k/2) = IN[  k/2 -> w-k/2) * K[0     -> k  )
 *      OUT[w-k/2 -> w    ) = IN[w-k/2 -> w    ) * K[0     -> k/2)
 * 
 * case: kernel length EVEN (ker. starts with 2 pix. out and ends with 1 pix. out)
 *      OUT[0       ->   k/2  ) = IN[0       ->   k/2  ) * K[  k/2 -> k    )
 *      OUT[  k/2   -> w-k/2+1) = IN[  k/2   -> w-k/2+1) * K[0     -> k    )
 *      OUT[w-k/2+1 -> w      ) = IN[w-k/2+1 -> w      ) * K[0     -> k/2+1)
 */
template<typename T>
void fft_h(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src, const Kernel1D<T>& kernel){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t new_channels = dest.getChannels();
    size_t ker_len = kernel.getLength();
    size_t k_l_2_a = (ker_len-1) / 2;
    size_t k_l_2_b = ker_len / 2;

    std::vector<T> priv_buf(new_channels, 0);

    // horizontal FT
    for(size_t i = 0; i < new_height; ++i){
        // part1
        for(size_t j = 0; j < k_l_2_a; ++j){

            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t k = k_l_2_a; k < ker_len; ++k){
                for(size_t c = 0; c < new_channels; ++c){
                    priv_buf[c] += src(i, j-k_l_2_a+k, c) * kernel(k, c);
                }
            }
            for(size_t c = 0; c < new_channels; ++c){
                dest(i, j, c) = static_cast<png_byte>(priv_buf[c]);
            }
        }
        // part2
        for(size_t j = k_l_2_a; j < new_width-k_l_2_b; ++j){
            
            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t k = 0; k < ker_len; ++k){
                for(size_t c = 0; c < new_channels; ++c){
                    priv_buf[c] += src(i, j-k_l_2_a+k, c) * kernel(k, c);
                }
            }
            for(size_t c = 0; c < new_channels; ++c){
                dest(i, j, c) = static_cast<png_byte>(priv_buf[c]);
            }
        }
        // part3
        for(size_t j = new_width-k_l_2_b; j < new_width; ++j){
            
            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t k = 0; k < k_l_2_a + (new_width-j); ++k){
                for(size_t c = 0; c < new_channels; ++c){
                    priv_buf[c] += src(i, j-k_l_2_a+k, c) * kernel(k, c);
                }
            }
            for(size_t c = 0; c < new_channels; ++c){
                dest(i, j, c) = static_cast<png_byte>(priv_buf[c]);
            }
        }
    }
}
/**
 * precomputes lookup table for kernel ranges
 * implementation for 'circular same fft'
 * for each index iterating over dest image,
 *      dest(n) += source(n-*(ker_inf_table++) : n+*(ker_sup_table++))
 * except for images we dont accumulate over the channels and keep
 * n_channels different values for each iteration step n
 * 
 * offset for X, indices for kernels
 * improve reduction: https://coderwall.com/p/gocbhg/openmp-improve-reduction-techniques
 */
template<typename T>
void fft_h2(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src, const Kernel1D<T>& kernel){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t new_channels = dest.getChannels();
    size_t ker_len = kernel.getLength();

    size_t k_l_2_a = (ker_len-1) / 2;
    size_t k_l_2_b = ker_len / 2;

    std::vector<size_t> offset_table(new_width, 0);
    std::vector<size_t> size_table(new_width, 0);
    for(size_t j = 0; j < k_l_2_a; ++j){
        offset_table[j] = j;
        size_table[j] = j + k_l_2_b + 1;
    }
    for(size_t j = k_l_2_a; j < new_width-k_l_2_b; ++j){
        offset_table[j] = k_l_2_a;
        size_table[j] = ker_len;
    }
    for(size_t j = new_width-k_l_2_b; j < new_width; ++j){
        offset_table[j] = k_l_2_a;
        size_table[j] = k_l_2_a + new_width - j;
    }

    std::vector<T> priv_buf(new_channels, 0);
    for(size_t i = 0; i < new_height; ++i){
        for(size_t j = 0; j < new_width; ++j){
            
            size_t l_offset = offset_table[j];
            size_t ker_size = size_table[j];
            std::fill(priv_buf.begin(), priv_buf.end(), 0);
            for(size_t k = 0; k < ker_size; ++k){
                for(size_t c = 0; c < new_channels; ++c){
                    priv_buf[c] += src(i, j-l_offset+k, c) * kernel(k_l_2_a-l_offset+k, c);
                }
            }
            for(size_t c = 0; c < new_channels; ++c){
                dest(i, j, c) = static_cast<png_byte>(priv_buf[c]);
            }
        }
    }
}

template<typename T>
void fft_h(BitMapRGBA& dest, const BitMapRGBA& src, const Kernel1D<T>& kernel){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t ker_size = kernel.getLength();
    size_t k_l_2_a = (ker_size-1) / 2;
    size_t k_l_2_b = ker_size / 2;

    // horizontal FT
    for(size_t i = 0; i < new_height; ++i){
        // part1
        for(size_t j = 0; j < k_l_2_a; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            T a = 0;
            for(size_t k = (k_l_2_a - j); k < ker_size; ++k){
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
                a += static_cast<T>((src(i, j-k_l_2_a+k) >> 24) & 0xff) * kernel(k, 3);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000) |
                         ((static_cast<uint32_t>(a) << 24) & 0xff000000);
        }
        // part2
        for(size_t j = k_l_2_a; j < new_width-k_l_2_b; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            T a = 0;
            for(size_t k = 0; k < ker_size; ++k){
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
                a += static_cast<T>((src(i, j-k_l_2_a+k) >> 24) & 0xff) * kernel(k, 3);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000) |
                         ((static_cast<uint32_t>(a) << 24) & 0xff000000);
        }
        // part3
        for(size_t j = new_width-k_l_2_b; j < new_width; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            T a = 0;
            for(size_t k = 0; k < k_l_2_a + (new_width-j); ++k){
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
                a += static_cast<T>((src(i, j-k_l_2_a+k) >> 24) & 0xff) * kernel(k, 3);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000) |
                         ((static_cast<uint32_t>(a) << 24) & 0xff000000);
        }
    }
}

template<typename T>
void fft_h2(BitMapRGBA& dest, const BitMapRGBA& src, const Kernel1D<T>& kernel){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();

    size_t ker_len = kernel.getLength();
    size_t k_l_2_a = (ker_len-1) / 2;

    ker_mask<true> km_w(ker_len, k_l_2_a, new_width);
    for(size_t i = 0; i < new_height; ++i){
        for(size_t j = 0; j < new_width; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            T a = 0;
            for(size_t k = 0; k < k_l_2_a; ++k){
                size_t offset = km_w(j, k);
                r += ( src(i, j-offset)        & 0xff) * kernel(k, 0);
                g += ((src(i, j-offset) >> 8 ) & 0xff) * kernel(k, 1);
                b += ((src(i, j-offset) >> 16) & 0xff) * kernel(k, 2);
                a += ((src(i, j-offset) >> 24) & 0xff) * kernel(k, 3);
            }
            for(size_t k = k_l_2_a; k < ker_len; ++k){
                size_t offset = km_w(j, k);
                r += ( src(i, j+offset)          & 0xff) * kernel(k, 0);
                g += ((src(i, j+offset) >> 8 )   & 0xff) * kernel(k, 1);
                b += ((src(i, j+offset) >> 16)   & 0xff) * kernel(k, 2);
                a += ((src(i, j+offset) >> 24)   & 0xff) * kernel(k, 3);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000) |
                         ((static_cast<uint32_t>(a) << 24) & 0xff000000);
        }
    }
}

template<typename T>
void fft_h(BitMapRGB& dest, const BitMapRGB& src, const Kernel1D<T>& kernel){
    
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();
    size_t ker_size = kernel.getLength();
    size_t k_l_2_a = (ker_size-1) / 2;
    size_t k_l_2_b = ker_size / 2;

    // horizontal FT
    for(size_t i = 0; i < new_height; ++i){
        // part1
        for(size_t j = 0; j < k_l_2_a; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            for(size_t k = (k_l_2_a - j); k < ker_size; ++k){
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000);
        }
        // part2
        for(size_t j = k_l_2_a; j < new_width-k_l_2_b; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            for(size_t k = 0; k < ker_size; ++k){
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000);
        }
        // part3
        for(size_t j = new_width-k_l_2_b; j < new_width; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            for(size_t k = 0; k < k_l_2_a + (new_width-j); ++k){
                //assert(j-k_l_2_a+k < new_width);
                r += static_cast<T>( src(i, j-k_l_2_a+k)        & 0xff) * kernel(k, 0);
                g += static_cast<T>((src(i, j-k_l_2_a+k) >> 8 ) & 0xff) * kernel(k, 1);
                b += static_cast<T>((src(i, j-k_l_2_a+k) >> 16) & 0xff) * kernel(k, 2);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000);
        }
    }
}

template<typename T>
void fft_h2(BitMapRGB& dest, const BitMapRGB& src, const Kernel1D<T>& kernel){
    size_t new_height = dest.getHeight();
    size_t new_width = dest.getWidth();

    size_t ker_len = kernel.getLength();
    size_t k_l_2_a = (ker_len-1) / 2;

    ker_mask<true> km_w(ker_len, k_l_2_a, new_width);
    for(size_t i = 0; i < new_height; ++i){
        for(size_t j = 0; j < new_width; ++j){
            T r = 0;
            T g = 0;
            T b = 0;
            for(size_t k = 0; k < k_l_2_a; ++k){
                size_t offset = km_w(j, k);
                r += ( src(i, j-offset)        & 0xff) * kernel(k, 0);
                g += ((src(i, j-offset) >> 8 ) & 0xff) * kernel(k, 1);
                b += ((src(i, j-offset) >> 16) & 0xff) * kernel(k, 2);
            }
            for(size_t k = k_l_2_a; k < ker_len; ++k){
                size_t offset = km_w(j, k);
                r += ( src(i, j+offset)          & 0xff) * kernel(k, 0);
                g += ((src(i, j+offset) >> 8 )   & 0xff) * kernel(k, 1);
                b += ((src(i, j+offset) >> 16)   & 0xff) * kernel(k, 2);
            }
            dest(i, j) = ( static_cast<uint32_t>(r)        & 0xff) |
                         ((static_cast<uint32_t>(g) << 8 ) & 0xff00) |
                         ((static_cast<uint32_t>(b) << 16) & 0xff0000);
        }
    }
}
/**
 * Fast Fourier Transform in 4 steps
 *      horizontal fft on input -> out
 *      transpose of out -> out2
 *      horizontal fft on out2 -> out3
 *      transpose of out3 -> out_final
*/
template<typename T>
void fft(PixelMap<png_byte>& dest, const PixelMap<png_byte>& src, const Kernel1D<T>& kernel){

    fft_h(dest, src, kernel);
    dest = dest.transpose();
    PixelMap<png_byte> temp = dest;
    fft_h(dest, temp, kernel);
    dest = dest.transpose();
}

template<typename T>
void fft(BitMapRGBA& dest, const BitMapRGBA& src, const Kernel1D<T>& kernel){

    fft_h(dest, src, kernel);
    dest = dest.transpose();
    BitMapRGBA temp = dest;
    fft_h(dest, temp, kernel);
    dest = dest.transpose();
}

template<typename T>
void fft(BitMapRGB& dest, const BitMapRGB& src, const Kernel1D<T>& kernel){
    
    fft_h(dest, src, kernel);
    dest = dest.transpose();
    BitMapRGB temp = dest;
    fft_h(dest, temp, kernel);
    dest = dest.transpose();
}

void DoG(PixelMap<png_byte>& dest,
        PixelMap<png_byte>& temp,
        const PixelMap<png_byte>& src,
        double sd1, double sd2){

    size_t k_w1 = std::ceil(sd1 * 6);
    size_t k_w2 = std::ceil(sd2 * 6);

    Kernel1D<float> ker1(k_w1, src.getChannels(), sd1);
    Kernel1D<float> ker2(k_w2, src.getChannels(), sd2);

    fft(temp, src, ker1);
    fft(dest, src, ker2);

    dest.subtract(temp);
}

void DoG(BitMapRGBA& dest,
        BitMapRGBA& temp,
        const BitMapRGBA& src,
        double sd1, double sd2){
    
    size_t k_w1 = std::ceil(sd1 * 6);
    size_t k_w2 = std::ceil(sd2 * 6);

    Kernel1D<float> ker1(k_w1, 4, sd1);
    Kernel1D<float> ker2(k_w2, 4, sd2);

    fft(temp, src, ker1);
    fft(dest, src, ker2);

    dest.subtract_rgb(temp);
}

void DoG(BitMapRGB& dest,
        BitMapRGB& temp,
        const BitMapRGB& src,
        double sd1, double sd2){
    
    size_t k_w1 = std::ceil(sd1 * 6);
    size_t k_w2 = std::ceil(sd2 * 6);

    Kernel1D<float> ker1(k_w1, 4, sd1);
    Kernel1D<float> ker2(k_w2, 4, sd2);

    fft(temp, src, ker1);
    fft(dest, src, ker2);

    dest.subtract(temp);
}