#pragma once

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <string_view>
#include <cstdarg>
#include <cassert>
#include <memory>
#include <limits.h>     // CHAR_BIT
// debug
#include <iostream>
// sudo apt-get install libpng-dev
// header path: /usr/include/libpng/png.h
#include <png.h>

#include "BitMap.h"

/**
 * http://libpng.org/pub/png/libpng-1.0.3-manual.html
 * http://www.libpng.org/pub/png/spec/1.2/PNG-Rationale.html#R.Byte-order
   Color    Allowed    Interpretation
   Type    Bit Depths
   
   0       1,2,4,8,16  Each pixel is a grayscale sample.
   
   2       8,16        Each pixel is an R,G,B triple.
   
   3       1,2,4,8     Each pixel is a palette index;
                       a PLTE chunk must appear.
   
   4       8,16        Each pixel is a grayscale sample,
                       followed by an alpha sample.
   
   6       8,16        Each pixel is an R,G,B triple,
                       followed by an alpha sample.
*/

typedef std::unique_ptr<FILE, int (*)(FILE*)> unique_fp;

static unique_fp make_unique_fp(const char* filename, const char* flags){
    return unique_fp(fopen(filename, flags), fclose);
}

struct PNG_utils{
    static constexpr int col_type[4] = {PNG_COLOR_TYPE_GRAY, 
                                        PNG_COLOR_TYPE_GRAY_ALPHA, 
                                        PNG_COLOR_TYPE_RGB, 
                                        PNG_COLOR_TYPE_RGB_ALPHA};
}; PNG_utils png_utils;

const png_dims getDims(std::string_view file_name){
    char header[8]; // max size that can be checked
    // open and test if png file
    auto fp = make_unique_fp(file_name.data(), "rb");
    if(!fp.get()) abort();

    if(!fread(header, 1, 8, fp.get())) abort();
        
    if(png_sig_cmp((png_const_bytep)header, 0, 8)) abort();
        
    // initializing
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr) abort();
        
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        abort();
    }
    png_infop end_info = png_create_info_struct(png_ptr);
    if(!end_info){
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
        abort();
    }
    if(setjmp(png_jmpbuf(png_ptr))){
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        abort();
    }
    png_init_io(png_ptr, fp.get());
    // lets libpng know there are some bytes missing (the 8 we read)
    png_set_sig_bytes(png_ptr, 8);

    // read all the file information up to the actual image data
    png_read_info(png_ptr, info_ptr);
    size_t height = png_get_image_height(png_ptr, info_ptr);
    size_t width = png_get_image_width(png_ptr, info_ptr);
    size_t channels = png_get_channels(png_ptr, info_ptr);
    size_t bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    png_byte color_type = png_get_color_type(png_ptr, info_ptr);

    png_read_update_info(png_ptr, info_ptr);

    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    png_ptr = NULL;
    info_ptr = NULL;
    end_info = NULL;
    return png_dims(height, width, channels, bit_depth, color_type);
}

void read_png_file(std::string_view file_name, PixelMap<png_byte>& pixelmap){
    char header[8]; // max size that can be checked
    // open and test if png file
    auto fp = make_unique_fp(file_name.data(), "rb");
    if(!fp.get()) abort();

    if(!fread(header, 1, 8, fp.get())) abort();
        
    if(png_sig_cmp((png_const_bytep)header, 0, 8)) abort();
        
    // initializing
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr) abort();
        
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        abort();
    }
    png_infop end_info = png_create_info_struct(png_ptr);
    if(!end_info){
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
        abort();
    }
    if(setjmp(png_jmpbuf(png_ptr))){
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        abort();
    }
    png_init_io(png_ptr, fp.get());
    // lets libpng know there are some bytes missing (the 8 we read)
    png_set_sig_bytes(png_ptr, 8);
        
    // read all the file information up to the actual image data
    png_read_info(png_ptr, info_ptr);
    size_t height = png_get_image_height(png_ptr, info_ptr);
    size_t width = png_get_image_width(png_ptr, info_ptr);
    size_t channels = png_get_channels(png_ptr, info_ptr);
    //size_t bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    //size_t rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    //png_byte color_type = png_get_color_type(png_ptr, info_ptr);
        
    //size_t n_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    // read file
    if(setjmp(png_jmpbuf(png_ptr))) abort();

    if(height != pixelmap.getHeight() || width != pixelmap.getWidth() || channels != pixelmap.getChannels())
        pixelmap = PixelMap<png_byte>(height, width, channels);        
    for(size_t y = 0; y < height; ++y)
        png_read_row(png_ptr, pixelmap.rowBegin(y), NULL);

    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    png_ptr = NULL;
    info_ptr = NULL;
    end_info = NULL;
}

template <typename COL_T>
void write_png_file(std::string_view file_name, const PixelMap<COL_T>& pixelmap){
    auto fp = make_unique_fp(file_name.data(), "wb");
    if(!fp.get()) abort();

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr) abort();
        
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        abort();
    }
    if(setjmp(png_jmpbuf(png_ptr))){
        png_destroy_write_struct(&png_ptr, &info_ptr);
        abort();
    }
    png_init_io(png_ptr, fp.get());
    /* write header */
    if(setjmp(png_jmpbuf(png_ptr))) abort();
    png_set_IHDR(png_ptr, info_ptr, 
                pixelmap.getWidth(), pixelmap.getHeight(),
                (CHAR_BIT*sizeof(COL_T)), png_utils.col_type[pixelmap.getChannels()-1], 
                PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
        
    // write bytes
    if(setjmp(png_jmpbuf(png_ptr))) abort();
    for(size_t y = 0; y < pixelmap.getHeight(); ++y)
        png_write_row(png_ptr, const_cast<png_const_bytep>(pixelmap.rowBegin(y)));
    // end write
    if(setjmp(png_jmpbuf(png_ptr))) abort();
    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);
    png_ptr = NULL;
    info_ptr = NULL;
}

void write_png_file(std::string_view file_name, const BitMapRGBA& bitmap){
    PixelMap<png_byte> pixelmap(bitmap);
    write_png_file(file_name, pixelmap);
}

void write_png_file(std::string_view file_name, const BitMapRGB& bitmap){
    PixelMap<png_byte> pixelmap(bitmap);
    write_png_file(file_name, pixelmap);
}