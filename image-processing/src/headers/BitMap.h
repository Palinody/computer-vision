#pragma once

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <string_view>
#include <cstdarg>
#include <cassert>
#include <memory>
// debug
#include <iostream>
#include <vector>

#include "PRNG.h"

enum DistrName{
    UNIFORM, 
    NORMAL 
};

struct png_dims{
    png_dims(size_t h, size_t w, size_t c, size_t bd, png_byte ct) : 
        height{ h }, width{ w }, channels{ c }, bit_depth{ bd },
        color_type(ct){}

    size_t height, width, channels, bit_depth;
    png_byte color_type;
};

// utility lookup maps structure to retrieve the correct shift/value for r-g-b-a pixels
struct RGBA_utils{
    static constexpr uint32_t rgba[4] = {0x000000ff, 
                                         0x0000ff00, 
                                         0x00ff0000, 
                                         0xff000000};
    static constexpr unsigned char shift[4] = {0, 8, 16, 24};
}; RGBA_utils rgba_utils;

#define CLIP(p, inf, sup) if(p < inf) { p = inf; } else if(p > sup) { p = sup; }

class BitMapRGBA;
class BitMapRGB;

template<typename T>
class PixelMap{
public:
    explicit
    PixelMap();
    PixelMap(size_t height, size_t width, size_t channels);
    PixelMap(png_dims&& dims);
    PixelMap(PixelMap<T>&& other) noexcept;
    PixelMap(const PixelMap<T>& other);
    PixelMap(const BitMapRGBA& bitmap);
    PixelMap(const BitMapRGB& bitmap);

    inline size_t getHeight() const;
    inline size_t getWidth() const;
    inline size_t getChannels() const;

    inline bool has_same_dims_as(const PixelMap<T>& other) const;
    PixelMap<T> transpose() const;
    PixelMap<T>& subtract(const PixelMap<T>& other);

    inline T& operator()(size_t i, size_t j, size_t c);
    inline const T& operator()(size_t i, size_t j, size_t c) const;
    
    inline T* begin();
    inline T* end();
    inline const T* begin() const;
    inline const T* end() const;

    inline T* rowBegin(size_t row);
    inline T* rowEnd(size_t row);
    inline const T* rowBegin(size_t row) const;
    inline const T* rowEnd(size_t row) const;

    inline T* pixBegin(size_t row, size_t col);
    inline T* pixEnd(size_t row, size_t col);
    inline const T* pixBegin(size_t row, size_t col) const;
    inline const T* pixEnd(size_t row, size_t col) const;

    PixelMap<T>& operator=(const PixelMap<T>& other);
    bool operator==(const PixelMap<T>& other) const;

    void print() const;

private:
    size_t _height = 0;
    size_t _width = 0;
    size_t _channels = 0;
    std::unique_ptr<T[]> _pixel_map;
};

class BitMapRGBA{
public:
    explicit
    BitMapRGBA();
    BitMapRGBA(size_t height, size_t width);
    BitMapRGBA(size_t height, size_t width, DistrName d_name);
    template<typename U, typename V>
    BitMapRGBA(size_t height, size_t width, DistrName d_name, U param1, V param2);
    BitMapRGBA(const BitMapRGBA& other);
    BitMapRGBA(BitMapRGBA&& other) noexcept;
    BitMapRGBA(const PixelMap<png_byte>& pixelmap);

    inline size_t getHeight() const;
    inline size_t getWidth() const;
    inline unsigned char 
    get(size_t i, size_t j, size_t c) const;

    template<typename rowIt>
    void copy_row(rowIt rowit, size_t row, size_t channels);
    BitMapRGBA transpose() const;
    BitMapRGBA& subtract_rgb(const BitMapRGBA& other);
    BitMapRGBA& subtract_rgba(const BitMapRGBA& other);

    inline uint32_t& operator()(size_t i, size_t j);
    inline const uint32_t& operator()(size_t i, size_t j) const;

    inline uint32_t* begin();
    inline uint32_t* end();
    inline const uint32_t* begin() const;
    inline const uint32_t* end() const;

    inline uint32_t* rowBegin(size_t row);
    inline uint32_t* rowEnd(size_t row);
    inline const uint32_t* rowBegin(size_t row) const;
    inline const uint32_t* rowEnd(size_t row) const;

    BitMapRGBA& operator=(const BitMapRGBA& other);
    //BitMapRGBA& operator=(const PixelMap<png_byte>& pixelmap);

    void print() const;
    void print_bitmap() const;
private:
    size_t _height = 0;
    size_t _width = 0;
    std::unique_ptr<uint32_t[]> _data;
};

class BitMapRGB{
public:
    explicit
    BitMapRGB();
    BitMapRGB(size_t height, size_t width);
    BitMapRGB(size_t height, size_t width, DistrName d_name);
    template<typename U, typename V>
    BitMapRGB(size_t height, size_t width, DistrName d_name, U param1, V param2);
    BitMapRGB(const BitMapRGB& other);
    BitMapRGB(BitMapRGB&& other) noexcept;
    BitMapRGB(const PixelMap<png_byte>& pixelmap);

    inline size_t getHeight() const;
    inline size_t getWidth() const;
    inline unsigned char 
    get(size_t i, size_t j, size_t c) const;

    template<typename rowIt>
    void copy_row(rowIt rowit, size_t row, size_t channels);
    BitMapRGB transpose() const;
    BitMapRGB& subtract(const BitMapRGB& other);

    inline uint32_t& operator()(size_t i, size_t j);
    inline const uint32_t& operator()(size_t i, size_t j) const;

    inline uint32_t* begin();
    inline uint32_t* end();
    inline const uint32_t* begin() const;
    inline const uint32_t* end() const;

    inline uint32_t* rowBegin(size_t row);
    inline uint32_t* rowEnd(size_t row);
    inline const uint32_t* rowBegin(size_t row) const;
    inline const uint32_t* rowEnd(size_t row) const;

    BitMapRGB& operator=(const BitMapRGB& other);
    //BitMapRGB& operator=(const PixelMap<png_byte>& pixelmap);
    //BitMapRGB& operator-=(const BitMapRGB other);

    void print() const;
    void print_bitmap() const;
private:
    size_t _height = 0;
    size_t _width = 0;
    std::unique_ptr<uint32_t[]> _data;
};

//////////////////////////////////////////////////////////////////////
//                                                                  //
//							 PIXELMAP IMPLEM                        //
//                                                                  //
//////////////////////////////////////////////////////////////////////

template<typename T>
PixelMap<T>::PixelMap() : 
    _height{ 0 }, _width{ 0 }, _channels{ 0 }{}

template<typename T>
PixelMap<T>::PixelMap(size_t height, size_t width, size_t channels) : 
    _height{ height }, _width{ width }, _channels{ channels },
    _pixel_map{ std::make_unique<T[]>(height*width*channels) }{}

template<typename T>
PixelMap<T>::PixelMap(png_dims&& dims) : 
        _height{ dims.height }, _width{ dims.width }, _channels{ dims.channels },
        _pixel_map{ std::make_unique<T[]>(dims.height*dims.width*dims.channels) }{}

template<typename T>
PixelMap<T>::PixelMap(PixelMap<T>&& other) noexcept:
        _height{ other._height }, _width{ other._width }, _channels{ other._channels },
        _pixel_map{ std::move(other._pixel_map) }{}

template<typename T>
PixelMap<T>::PixelMap(const PixelMap<T>& other){
    *this = other;
}

template<typename T>
PixelMap<T>::PixelMap(const BitMapRGBA& bitmap) :
    _height{ bitmap.getHeight() }, _width{ bitmap.getWidth() }, _channels{ 4 },
    _pixel_map{ std::make_unique<T[]>(bitmap.getHeight()*bitmap.getWidth()*4) }{
        
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            (*this)(i, j, 0) = ( bitmap(i, j)       & 0xff);
            (*this)(i, j, 1) = ((bitmap(i, j) >> 8 )& 0xff);
            (*this)(i, j, 2) = ((bitmap(i, j) >> 16)& 0xff);
            (*this)(i, j, 3) = ((bitmap(i, j) >> 24)& 0xff);
        }
    }
}

template<typename T>
PixelMap<T>::PixelMap(const BitMapRGB& bitmap) :
    _height{ bitmap.getHeight() }, _width{ bitmap.getWidth() }, _channels{ 3 },
    _pixel_map{ std::make_unique<T[]>(bitmap.getHeight()*bitmap.getWidth()*3) }{
        
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            (*this)(i, j, 0) = ( bitmap(i, j)       & 0xff);
            (*this)(i, j, 1) = ((bitmap(i, j) >> 8 )& 0xff);
            (*this)(i, j, 2) = ((bitmap(i, j) >> 16)& 0xff);
        }
    }
}

template<typename T>
size_t PixelMap<T>::getHeight() const{ return _height; }

template<typename T>
size_t PixelMap<T>::getWidth() const{ return _width; }

template<typename T>
size_t PixelMap<T>::getChannels() const{ return _channels; }

template<typename T>
T& PixelMap<T>::operator()(size_t i, size_t j, size_t c){ return _pixel_map[c+(j+i*_width)*_channels]; }
template<typename T>
const T& PixelMap<T>::operator()(size_t i, size_t j, size_t c) const{ return _pixel_map[c+(j+i*_width)*_channels]; }

template<typename T>
T* PixelMap<T>::begin() { return _pixel_map.get(); }
template<typename T>
T* PixelMap<T>::end() { return begin() + _width*_height*_channels; }

template<typename T>
const T* PixelMap<T>::begin() const { return _pixel_map.get(); }
template<typename T>
const T* PixelMap<T>::end() const { return begin() + _width*_height*_channels; }

template<typename T>
T* PixelMap<T>::rowBegin(size_t row) { return begin() + row*_width*_channels; }
template<typename T>
T* PixelMap<T>::rowEnd(size_t row) { return rowBegin(row) + _width*_channels; }
template<typename T>
const T* PixelMap<T>::rowBegin(size_t row) const { return begin() + row*_width*_channels; }
template<typename T>
const T* PixelMap<T>::rowEnd(size_t row) const { return rowBegin(row) + _width*_channels; }

template<typename T>
T* PixelMap<T>::pixBegin(size_t row, size_t col) { return rowBegin(row) + col*_channels; }
template<typename T>
T* PixelMap<T>::pixEnd(size_t row, size_t col) { return pixBegin(row, col)  + _channels; }
template<typename T>
const T* PixelMap<T>::pixBegin(size_t row, size_t col) const { return rowBegin(row) + col*_channels; }
template<typename T>
const T* PixelMap<T>::pixEnd(size_t row, size_t col) const { return pixBegin(row, col)  + _channels; }

template<typename T>
PixelMap<T>& PixelMap<T>::operator=(const PixelMap<T>& other){
    if(_height != other.getHeight() || _width != other.getWidth() || _channels != other.getChannels()){
        _height = other.getHeight();
        _width = other.getWidth();
        _channels = other.getChannels();
        _pixel_map.reset(nullptr);
        _pixel_map = std::make_unique<T[]>(_height*_width*_channels);
    }
    T* it_dest = begin();
    const T* it_src = other.begin();
    for(; it_dest != end(); ++it_dest, ++it_src) *it_dest = *it_src;

    return *this;
}

template<typename T>
bool PixelMap<T>::operator==(const PixelMap<T>& other) const{
    if(has_same_dims_as(other))
        return false;
    return std::equal(begin(), other.begin(), other.end());
}

template<typename T>
bool PixelMap<T>::has_same_dims_as(const PixelMap<T>& other) const{
    if(_height != other.getHeight() || _width != other.getWidth() || _channels != other.getChannels())
        return false;
    return true;
}

template<typename T>
PixelMap<T> PixelMap<T>::transpose() const{
    
    PixelMap<T> res(_width, _height, _channels);

    for(size_t y = 0; y < _height; ++y){
        for(size_t x = 0; x < _width; ++x){
            for(size_t c = 0; c < _channels; ++c){
                res(x, y, c) = (*this)(y, x, c);
            }
        }
    }
    return res;
}
// if min val negative shift all values 
/**
 * normalization:
 * x_new = (x - min)*(max_des - min_des)/(max-min)+min_des
*/
template<typename T>
PixelMap<T>& PixelMap<T>::subtract(const PixelMap<T>& other){

    assert(has_same_dims_as(other));

    std::vector<uint16_t> min_val(_channels);
    std::vector<uint16_t> max_val(_channels);
    #ifdef _OPENMP
        #pragma omp for simd
    #endif
    for(size_t c = 0; c < _channels; ++c) min_val[c] = 510; // 2*255
    #ifdef _OPENMP
        #pragma omp for simd
    #endif
    for(size_t c = 0; c < _channels; ++c) max_val[c] = 0;

    std::vector<uint16_t> buf(_width*_height*_channels);
    #ifdef _OPENMP
        #pragma omp parallel for simd collapse(3)
    #endif
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            for(size_t c = 0; c < _channels; ++c){
                buf[c+(j+i*_width)*_channels] = (255 + (*this)(i, j, c)) - static_cast<uint32_t>(other(i, j, c));
                if(buf[c+(j+i*_width)*_channels] < min_val[c]) min_val[c] = buf[c+(j+i*_width)*_channels];
                else
                if(buf[c+(j+i*_width)*_channels] > max_val[c]) max_val[c] = buf[c+(j+i*_width)*_channels];
            }
        }
    }
    #ifdef _OPENMP
        #pragma omp parallel for simd collapse(3)
    #endif
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            for(size_t c = 0; c < _channels; ++c){
                float fact = 255/(max_val[c]-min_val[c]+1e-10f/* == 0 ? 1e-10f : max_val[c]-min_val[c]*/);
                (*this)(i, j, c) = static_cast<png_byte>((buf[c+(j+i*_width)*_channels] - min_val[c]) * fact);
            }
        }
    }
    return *this;
}

template<typename T>
void PixelMap<T>::print() const{
    const T* it = begin();
    for(size_t y = 0; y < _height; ++y){
        for(size_t x = 0; x < _width; ++x){
            //const T* it = pixBegin(y, x);
            std::cout << "{";
            for(size_t c = 0; c < _channels; ++c){
                std::cout << *static_cast<const T*>(it++) << ", ";
            }
            std::cout << "}, ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////
//                                                                  //
//							 BITMAP IMPLEM (RGBA)                   //
//                                                                  //
//////////////////////////////////////////////////////////////////////

BitMapRGBA::BitMapRGBA() : 
    _height{ 0 }, _width{ 0 }{}

BitMapRGBA::BitMapRGBA(size_t height, size_t width) : 
    _height{ height }, _width{ width },
    _data{ std::make_unique<uint32_t[]>(height*width) }{}

/**
 * Constructor with default arguments that may have different types
 * d_name == UNIFORM -> pixel in [0, 255]
 * d_name == NORMAL -> n_trials = height*width*channels | probability = 0.5
*/
BitMapRGBA::BitMapRGBA(size_t height, size_t width, DistrName d_name) : 
    BitMapRGBA(height, width, d_name, ((d_name==NORMAL) ? 255 : 0), ((d_name==NORMAL) ? 0.5f : 255)){}


template<typename U, typename V>
BitMapRGBA::BitMapRGBA(size_t height, size_t width, DistrName d_name, U param1, V param2) : 
    _height{ height }, _width{ width },
    _data{ std::make_unique<uint32_t[]>(height*width) }{ 

    std::unique_ptr<distribution<uint32_t>> distr;
    switch(d_name){
        case UNIFORM:{
            distr = std::make_unique<uniform_dist<uint32_t>>(param1, param2);
            break;
        }
        case NORMAL:{
            distr = std::make_unique<normal_dist<uint32_t>>(param1, param2);
            break;
        }
        default:{
            printf("default bitmap noise: uniform\n");
            distr = std::make_unique<uniform_dist<uint32_t>>(param1, param2);
            break;
        }
    };

    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            _data[j+i*_width] = (((*distr)() << 0 )& 0x000000ff) |
                                (((*distr)() << 8 )& 0x0000ff00) |
                                (((*distr)() << 16)& 0x00ff0000) |
                                (((*distr)() << 24)& 0xff000000);
        }
    }
}

BitMapRGBA::BitMapRGBA(const BitMapRGBA& other){
    *this = other;
}

BitMapRGBA::BitMapRGBA(BitMapRGBA&& other) noexcept:
    _height{ other._height }, _width{ other._width },
    _data{ std::move(other._data) }{}

BitMapRGBA::BitMapRGBA(const PixelMap<png_byte>& pixelmap) :
    BitMapRGBA(pixelmap.getHeight(), pixelmap.getWidth()){
    
    size_t channels = pixelmap.getChannels();
    switch(channels){
        case 1:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j);
                    *(it_bm++) = (                                   0xff00) |
                                 (static_cast<uint32_t>(*pixel_pm) & 0x00ff);
                }
            }
            break;
        }
        case 2:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j);
                    *(it_bm++) = (static_cast<uint32_t>(*(pixel_pm+1) << 8) & 0x0000ff00) |
                                 (static_cast<uint32_t>(* pixel_pm    << 0) & 0x000000ff);
                }
            }
            break;
        }
        case 3:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j);
                    *(it_bm++) = (                                               0xff000000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+2)) << 16) & 0x00ff0000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+1)) << 8 ) & 0x0000ff00) |
                                 ( static_cast<uint32_t>(* pixel_pm   )        & 0x000000ff);
                }
            }
            break;
        }
        case 4:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j);
                    *(it_bm++) = ((static_cast<uint32_t>(*(pixel_pm+3)) << 24) & 0xff000000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+2)) << 16) & 0x00ff0000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+1)) << 8 ) & 0x0000ff00) |
                                 ( static_cast<uint32_t>(* pixel_pm   )        & 0x000000ff);
                }
            }
            break;
        }
        default:{
            printf("Must provide (1-4) channels for Bitmap construction.\n");
            abort();
        }
    }
}

size_t BitMapRGBA::getHeight() const{ return _height; }
size_t BitMapRGBA::getWidth() const{ return _width; }
unsigned char 
BitMapRGBA::get(size_t i, size_t j, size_t c) const{ return (_data[j+i*_width] >> (8*c) & 0xff); }

uint32_t& BitMapRGBA::operator()(size_t i, size_t j){ return _data[j+i*_width]; }
const uint32_t& BitMapRGBA::operator()(size_t i, size_t j) const{ return _data[j+i*_width]; }

uint32_t* BitMapRGBA::begin(){ return _data.get(); }
uint32_t* BitMapRGBA::end(){ return begin() + _height*_width; }
const uint32_t* BitMapRGBA::begin() const{ return _data.get(); }
const uint32_t* BitMapRGBA::end() const{ return begin() + _height*_width; }

uint32_t* BitMapRGBA::rowBegin(size_t row){ return begin() + row*_height*_width; }
uint32_t* BitMapRGBA::rowEnd(size_t row){ return rowBegin(row) + _width; }
const uint32_t* BitMapRGBA::rowBegin(size_t row) const{ return begin() + row*_height*_width; }
const uint32_t* BitMapRGBA::rowEnd(size_t row) const{ return rowBegin(row) + _width; }

template<typename rowIt>
void BitMapRGBA::copy_row(rowIt rowit, size_t row, size_t channels){
    
    //assert(row < _height && (*rowit)->getWidth() == _width);
    for(size_t j = 0; j < _width; ++j){
        for(size_t c = 0; c < channels; ++c){
            _data[j+row*_width] |= (*(rowit+(c+j*channels)) << rgba_utils.shift[c] & rgba_utils.rgba[c]);
        }
    }
}

BitMapRGBA BitMapRGBA::transpose() const{

    BitMapRGBA res(_width, _height);
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            res(j, i) = _data[j + i*_width];
        }
    }
    return res;
}
/*
 * shift all values by min_value to avoid negative values 
 *
 * normalization:
 * x_new = (x - min)*(max_des - min_des)/(max-min)+min_des
 *
 * each rgba max value "spreads" in uint64_t integer
 * max possible value: 2 * 255 = 510 = 0x01fe
*/
BitMapRGBA& BitMapRGBA::subtract_rgb(const BitMapRGBA& other){

    uint64_t min_val = 0x000001fe01fe01fe;
    uint64_t max_val = 0x0;
    std::vector<uint64_t> buf(_width*_height, 0);
    for(size_t n = 0; n < _width*_height; ++n) buf[n] = 0x000000ff00ff00ff;
    
    const uint32_t *this_it = begin();
    const uint32_t *other_it = other.begin();
    for(size_t n = 0; n < _width*_height; ++this_it, ++other_it, ++n){
        
        buf[n] += ((static_cast<uint64_t>(*this_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*this_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*this_it & 0xff));
        
        buf[n] -= ((static_cast<uint64_t>(*other_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*other_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*other_it & 0xff));

        min_val = ((buf[n] & 0x01fe00000000) < (min_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (min_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) < (min_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (min_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) < (min_val & 0x01fe) ? (buf[n] & 0x01fe) : (min_val & 0x01fe));
        
        max_val = ((buf[n] & 0x01fe00000000) > (max_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (max_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) > (max_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (max_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) > (max_val & 0x01fe) ? (buf[n] & 0x01fe) : (max_val & 0x01fe));
    }
    uint64_t max_min = max_val - min_val;
    // norm. factor for each col. channel (diff. can still be 0x1fe at this point)
    float b_fact = 255.0f / ((max_min & 0x1fe00000000) >> 32);
    float g_fact = 255.0f / ((max_min & 0x1fe0000) >> 16);
    float r_fact = 255.0f / ( max_min & 0x1fe);
    
    uint32_t *this_it_update = begin();
    for(size_t n = 0; n < _width*_height; ++n, ++this_it_update){
        //res = (buf[n] - min_val) * fact;
        buf[n] -= min_val;
        // shifting the bits so that they line up in proper uint32_t position. Alpha channel is copied
        *this_it_update = ((*this_it_update) & 0xff000000) |
                          ((static_cast<uint32_t>(((buf[n] & 0xff00000000) >> 32) * b_fact)/* & 0xff*/) << 16) |
                          ((static_cast<uint32_t>(((buf[n] & 0x0000ff0000) >> 16) * g_fact)/* & 0xff*/) << 8 ) |
                          ((static_cast<uint32_t>(((buf[n] & 0x00000000ff)      ) * r_fact)/* & 0xff*/)      );
    }
    return *this;
}

BitMapRGBA& BitMapRGBA::subtract_rgba(const BitMapRGBA& other){

    uint64_t min_val = 0x01fe01fe01fe01fe;
    uint64_t max_val = 0x0;
    std::vector<uint64_t> buf(_width*_height, 0);
    for(size_t n = 0; n < _width*_height; ++n) buf[n] = 0x00ff00ff00ff00ff;
    
    const uint32_t *this_it = begin();
    const uint32_t *other_it = other.begin();
    for(size_t n = 0; n < _width*_height; ++n, ++this_it, ++other_it){
        
        buf[n] += ((static_cast<uint64_t>(*this_it & 0xff000000)) << 24) |
                  ((static_cast<uint64_t>(*this_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*this_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*this_it & 0xff));
        
        buf[n] -= ((static_cast<uint64_t>(*other_it & 0xff000000)) << 24) |
                  ((static_cast<uint64_t>(*other_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*other_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*other_it & 0xff));

        min_val = ((buf[n] & 0x01fe000000000000) < (min_val & 0x01fe000000000000) ? (buf[n] & 0x01fe000000000000) : (min_val & 0x01fe000000000000)) |
                  ((buf[n] & 0x01fe00000000) < (min_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (min_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) < (min_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (min_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) < (min_val & 0x01fe) ? (buf[n] & 0x01fe) : (min_val & 0x01fe));
        
        max_val = ((buf[n] & 0x01fe000000000000) > (max_val & 0x01fe000000000000) ? (buf[n] & 0x01fe000000000000) : (max_val & 0x01fe000000000000)) |
                  ((buf[n] & 0x01fe00000000) > (max_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (max_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) > (max_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (max_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) > (max_val & 0x01fe) ? (buf[n] & 0x01fe) : (max_val & 0x01fe));
    }
    uint64_t max_min = max_val - min_val;
    // norm. factor for each col. channel (diff. can still be 0x1fe at this point)
    float a_fact = 255.0f / ((max_min & 0x1fe000000000000) >> 48);
    float b_fact = 255.0f / ((max_min & 0x1fe00000000) >> 32);
    float g_fact = 255.0f / ((max_min & 0x1fe0000) >> 16);
    float r_fact = 255.0f / ( max_min & 0x1fe);
    
    uint32_t *this_it_update = begin();
    for(size_t n = 0; n < _width*_height; ++n, ++this_it_update){
        //res = (buf[n] - min_val) * fact;
        buf[n] -= min_val;
        // shifting the bits so that they line up in proper uint32_t position
        *this_it_update = ((static_cast<uint32_t>(((buf[n] & 0xff0000000000) >> 48) * a_fact)/* & 0xff*/) << 24) |
                          ((static_cast<uint32_t>(((buf[n] & 0x00ff00000000) >> 32) * b_fact)/* & 0xff*/) << 16) |
                          ((static_cast<uint32_t>(((buf[n] & 0x000000ff0000) >> 16) * g_fact)/* & 0xff*/) << 8 ) |
                          ((static_cast<uint32_t>(((buf[n] & 0x0000000000ff)      ) * r_fact)/* & 0xff*/)      );
    }
    return *this;
}

BitMapRGBA& BitMapRGBA::operator=(const BitMapRGBA& other){
    if(_height != other.getHeight() || _width != other.getWidth()){
        _height = other.getHeight();
        _width = other.getWidth();
        _data.reset(nullptr);
        _data = std::make_unique<uint32_t[]>(_height*_width);
    }
    uint32_t* it_dest = begin();
    const uint32_t* it_src = other.begin();
    for(; it_dest != end(); ++it_dest, ++it_src) *it_dest = *it_src;

    return *this;
}
/*
BitMapRGBA& BitMapRGBA::operator=(const PixelMap<png_byte>& pixelmap){
    
    *this = BitMapRGBA(pixelmap);
    return *this;
}
*/
void BitMapRGBA::print() const{
    
    printf("red = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 0 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 0 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
    printf("green = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 8 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 8 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
    printf("blue = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 16 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 16 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
    printf("alpha = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 24 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 24 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
}

void BitMapRGBA::print_bitmap() const{
    
    printf("red = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)((_data[j+i*_width] >> 0) & 0xff));    
            else
                printf("%d,", (int)((_data[j+i*_width] >> 0) & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    printf("green = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)((_data[j+i*_width] >> 8) & 0xff));    
            else
                printf("%d,", (int)((_data[j+i*_width] >> 8) & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    printf("blue = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)((_data[j+i*_width] >> 16) & 0xff));    
            else
                printf("%d,", (int)((_data[j+i*_width] >> 16) & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    /*
    printf("alpha = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)((_data[j+i*_width] >> 24) & 0xff));    
            else
                printf("%d,", (int)((_data[j+i*_width] >> 24) & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    */
}

//////////////////////////////////////////////////////////////////////
//                                                                  //
//							 BITMAP IMPLEM (RGB)                   //
//                                                                  //
//////////////////////////////////////////////////////////////////////

BitMapRGB::BitMapRGB() : 
    _height{ 0 }, _width{ 0 }{}

BitMapRGB::BitMapRGB(size_t height, size_t width) : 
    _height{ height }, _width{ width },
    _data{ std::make_unique<uint32_t[]>(height*width) }{}

/**
 * Constructor with default arguments that may have different types
 * d_name == UNIFORM -> pixel in [0, 255]
 * d_name == NORMAL -> n_trials = height*width*channels | probability = 0.5
*/
BitMapRGB::BitMapRGB(size_t height, size_t width, DistrName d_name) : 
    BitMapRGB(height, width, d_name, ((d_name==NORMAL) ? 255 : 0), ((d_name==NORMAL) ? 0.5f : 255)){}


template<typename U, typename V>
BitMapRGB::BitMapRGB(size_t height, size_t width, DistrName d_name, U param1, V param2) : 
    _height{ height }, _width{ width },
    _data{ std::make_unique<uint32_t[]>(height*width) }{ 

    std::unique_ptr<distribution<uint32_t>> distr;
    switch(d_name){
        case UNIFORM:{
            distr = std::make_unique<uniform_dist<uint32_t>>(param1, param2);
            break;
        }
        case NORMAL:{
            distr = std::make_unique<normal_dist<uint32_t>>(param1, param2);
            break;
        }
        default:{
            printf("default bitmap noise: uniform\n");
            distr = std::make_unique<uniform_dist<uint32_t>>(param1, param2);
            break;
        }
    };

    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            _data[j+i*_width] = (( *distr)()       &     0xff) |
                                (((*distr)() << 8 )&   0xff00) |
                                (((*distr)() << 16)& 0xff0000);
        }
    }
}

BitMapRGB::BitMapRGB(const BitMapRGB& other){
    *this = other;
}

BitMapRGB::BitMapRGB(BitMapRGB&& other) noexcept:
    _height{ other._height }, _width{ other._width },
    _data{ std::move(other._data) }{}

BitMapRGB::BitMapRGB(const PixelMap<png_byte>& pixelmap) :
    BitMapRGB(pixelmap.getHeight(), pixelmap.getWidth()){
    
    size_t channels = pixelmap.getChannels();
    switch(channels){
        case 1:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j); // first channel
                    *(it_bm++) = (static_cast<uint32_t>(*pixel_pm) & 0xff);
                                 
                }
            }
            break;
        }
        case 2:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j); // first and alpha channels
                    *(it_bm++) = (static_cast<uint32_t>(*(pixel_pm+1) <<  8) & 0xff00) |
                                 (static_cast<uint32_t>(* pixel_pm         ) & 0x00ff);
                }
            }
            break;
        }
        case 3:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j); // rgb channels
                    *(it_bm++) = ((static_cast<uint32_t>(*(pixel_pm+2)) << 16) & 0xff0000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+1)) << 8 ) & 0x00ff00) |
                                 ( static_cast<uint32_t>(* pixel_pm   )        & 0x0000ff);
                                 
                }
            }
            break;
        }
        case 4:{
            uint32_t *it_bm = begin();
            for(size_t i = 0; i < _height; ++i){
                for(size_t j = 0; j < _width; ++j){
                    const png_byte *pixel_pm = pixelmap.pixBegin(i, j); // skip alpha channel
                    *(it_bm++) = ((static_cast<uint32_t>(*(pixel_pm+2)) << 16) & 0xff0000) |
                                 ((static_cast<uint32_t>(*(pixel_pm+1)) << 8 ) & 0x00ff00) |
                                 ( static_cast<uint32_t>(* pixel_pm   )        & 0x0000ff);
                }
            }
            break;
        }
        default:{
            printf("Must provide (1-4) channels for Bitmap construction.\n");
            abort();
        }
    }
}

size_t BitMapRGB::getHeight() const{ return _height; }
size_t BitMapRGB::getWidth() const{ return _width; }
unsigned char 
BitMapRGB::get(size_t i, size_t j, size_t c) const{ return (_data[j+i*_width] >> (8*c) & 0xff); }

uint32_t& BitMapRGB::operator()(size_t i, size_t j){ return _data[j+i*_width]; }
const uint32_t& BitMapRGB::operator()(size_t i, size_t j) const{ return _data[j+i*_width]; }

uint32_t* BitMapRGB::begin(){ return _data.get(); }
uint32_t* BitMapRGB::end(){ return begin() + _height*_width; }
const uint32_t* BitMapRGB::begin() const{ return _data.get(); }
const uint32_t* BitMapRGB::end() const{ return begin() + _height*_width; }

uint32_t* BitMapRGB::rowBegin(size_t row){ return begin() + row*_height*_width; }
uint32_t* BitMapRGB::rowEnd(size_t row){ return rowBegin(row) + _width; }
const uint32_t* BitMapRGB::rowBegin(size_t row) const{ return begin() + row*_height*_width; }
const uint32_t* BitMapRGB::rowEnd(size_t row) const{ return rowBegin(row) + _width; }

BitMapRGB& BitMapRGB::operator=(const BitMapRGB& other){
    if(_height != other.getHeight() || _width != other.getWidth()){
        _height = other.getHeight();
        _width = other.getWidth();
        _data.reset(nullptr);
        _data = std::make_unique<uint32_t[]>(_height*_width);
    }
    uint32_t* it_dest = begin();
    const uint32_t* it_src = other.begin();
    for(; it_dest != end(); ++it_dest, ++it_src) *it_dest = *it_src;

    return *this;
}
/*
BitMapRGB& BitMapRGB::operator=(const PixelMap<png_byte>& pixelmap){
    
    *this = BitMapRGB(pixelmap);
    return *this;
}
*/
template<typename rowIt>
void BitMapRGB::copy_row(rowIt rowit, size_t row, size_t channels){
    assert(channels < 4);
    //assert(row < _height && (*rowit)->getWidth() == _width);
    for(size_t j = 0; j < _width; ++j){
        for(size_t c = 0; c < 3; ++c){
            _data[j+row*_width] |= (*(rowit+(c+j*channels)) << rgba_utils.shift[c] & rgba_utils.rgba[c]);
        }
    }
}

BitMapRGB BitMapRGB::transpose() const{

    BitMapRGB res(_width, _height);
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            res(j, i) = _data[j + i*_width];
        }
    }
    return res;
}

BitMapRGB& BitMapRGB::subtract(const BitMapRGB& other){

    uint64_t min_val = 0x01fe01fe01fe;
    uint64_t max_val = 0x0;
    std::vector<uint64_t> buf(_width*_height, 0);
    for(size_t n = 0; n < _width*_height; ++n) buf[n] = 0x00ff00ff00ff;
    
    const uint32_t *this_it = begin();
    const uint32_t *other_it = other.begin();
    for(size_t n = 0; n < _width*_height; ++this_it, ++other_it, ++n){
        
        buf[n] += ((static_cast<uint64_t>(*this_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*this_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*this_it & 0xff));
        
        buf[n] -= ((static_cast<uint64_t>(*other_it & 0xff0000)) << 16) |
                  ((static_cast<uint64_t>(*other_it & 0xff00)) << 8) |
                  ( static_cast<uint64_t>(*other_it & 0xff));

        min_val = ((buf[n] & 0x01fe00000000) < (min_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (min_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) < (min_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (min_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) < (min_val & 0x01fe) ? (buf[n] & 0x01fe) : (min_val & 0x01fe));
        
        max_val = ((buf[n] & 0x01fe00000000) > (max_val & 0x01fe00000000) ? (buf[n] & 0x01fe00000000) : (max_val & 0x01fe00000000)) |
                  ((buf[n] & 0x01fe0000) > (max_val & 0x01fe0000) ? (buf[n] & 0x01fe0000) : (max_val & 0x01fe0000)) |
                  ((buf[n] & 0x01fe) > (max_val & 0x01fe) ? (buf[n] & 0x01fe) : (max_val & 0x01fe));
    }
    uint64_t max_min = max_val - min_val;
    // norm. factor for each col. channel (diff. can still be 0x1fe at this point)
    float b_fact = 255.0f / ((max_min & 0x1fe00000000) >> 32);
    float g_fact = 255.0f / ((max_min & 0x1fe0000) >> 16);
    float r_fact = 255.0f / ( max_min & 0x1fe);
    
    uint32_t *this_it_update = begin();
    for(size_t n = 0; n < _width*_height; ++n, ++this_it_update){
        //res = (buf[n] - min_val) * fact;
        buf[n] -= min_val;
        // shifting the bits so that they line up in proper uint32_t position. Alpha channel is copied
        *this_it_update = ((static_cast<uint32_t>(((buf[n] & 0xff00000000) >> 32) * b_fact)/* & 0xff*/) << 16) |
                          ((static_cast<uint32_t>(((buf[n] & 0x0000ff0000) >> 16) * g_fact)/* & 0xff*/) << 8 ) |
                          ((static_cast<uint32_t>(((buf[n] & 0x00000000ff)      ) * r_fact)/* & 0xff*/)      );
    }
    return *this;
}

void BitMapRGB::print() const{
    
    printf("red = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 0 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 0 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
    printf("green = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 8 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 8 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
    printf("blue = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d]", (int)(_data[j+i*_width] >> 16 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 16 & 0xff));
        }
        //printf("\n");
    }
    printf("\n");
}

void BitMapRGB::print_bitmap() const{
    
    printf("red = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)(_data[j+i*_width] >> 0 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 0 & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    printf("green = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)(_data[j+i*_width] >> 8 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 8 & 0xff));
        }
        printf("\n");
    }
    printf("\n");
    printf("blue = [");
    for(size_t i = 0; i < _height; ++i){
        for(size_t j = 0; j < _width; ++j){
            if(i == _height-1 && j == _width-1)
                printf("%d],", (int)(_data[j+i*_width] >> 16 & 0xff));    
            else
                printf("%d,", (int)(_data[j+i*_width] >> 16 & 0xff));
        }
        printf("\n");
    }
    printf("\n");
}