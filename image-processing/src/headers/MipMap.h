#pragma once

#include "ImageProcessing.h"
#include <deque>
#include <type_traits>
#include <algorithm>
#include <vector>

/**
 * T: image type
 * P: element of channel type (currently no effect if bitmap)
*/
//template<typename T, typename P = png_byte>
template<typename T>
class MipMap{
public:
    using imgType = std::conditional_t<
    std::is_same_v<T, PixelMap<png_byte>>, 
        PixelMap<png_byte>, 
        std::conditional_t<
            std::is_same_v<T, BitMapRGB>, 
                BitMapRGB, 
                BitMapRGBA
        >
    >;
    /** 
     * is_same_v used because PixelMap and BitMaps take 
     * a diff. numb. of arguments to construct an object
    */
    // used for PixelMap<png_byte>
    MipMap(const T& img) requires(std::is_same_v<T, PixelMap<png_byte>>);
    // used for all the BitMaps
    MipMap(const T& img) requires(!std::is_same_v<T, PixelMap<png_byte>>);

    inline size_t getOctaves() const;
    void get_interpolated_img(float ratio) const requires(std::is_same_v<T, PixelMap<png_byte>>);
    void get_interpolated_img(float ratio) const requires(!std::is_same_v<T, PixelMap<png_byte>>);

    inline T& operator[](size_t idx);
    inline const T& operator[](size_t idx) const;

    void savefig(std::string_view filename) const requires(std::is_same_v<T, PixelMap<png_byte>>);
    void savefig(std::string_view filename) const requires(!std::is_same_v<T, PixelMap<png_byte>>);

public:
    std::deque<T> _mipmap;
    size_t _octaves = 0;
};

template<typename T>
MipMap<T>::MipMap(const T& img) requires(std::is_same_v<T, PixelMap<png_byte>>) 
        : _mipmap{ img }, _octaves{ 1 } {
    
    size_t shortest_dim = std::min(img.getHeight(), img.getWidth());
    assert(shortest_dim > 1);
    size_t n_img = 1 + static_cast<size_t>(std::log2(static_cast<float>(shortest_dim)) / std::log2(2.0f));
    
    _mipmap.resize(n_img);
    imgType res = img;
    std::deque<size_t> dims_h = { _mipmap.at(0).getHeight() };
    std::deque<size_t> dims_w = { _mipmap.at(0).getWidth() };
    do{
        imgType temp(res.getHeight()/2, res.getWidth()/2, res.getChannels());
        rescale::bicubic(temp, res);
        _mipmap.emplace_back(temp);
        res = temp;

        dims_h.emplace_back(res.getHeight());
        dims_w.emplace_back(res.getWidth());
        ++_octaves;
    }while(res.getHeight() > 1 && res.getWidth() > 1);

    printf("stored height dimensions\n");
    for(const auto& h : dims_h) printf("%ld, ", h); printf("\n");
    for(const auto& w : dims_w) printf("%ld, ", w); printf("\n");
}

template<typename T>
MipMap<T>::MipMap(const T& img) requires(!std::is_same_v<T, PixelMap<png_byte>>) 
        : _mipmap{ img }, _octaves{ 1 } {
    
    size_t shortest_dim = std::min(img.getHeight(), img.getWidth());
    assert(shortest_dim > 1);
    size_t n_img = 1 + static_cast<size_t>(std::log2(static_cast<float>(shortest_dim)) / std::log2(2.0f));

    _mipmap.resize(n_img);
    imgType res = img;
    do{
        imgType temp(res.getHeight()/2, res.getWidth()/2);
        rescale::bicubic(temp, res);
        _mipmap.emplace_back(temp);
        res = temp;
        ++_octaves;
    }while(res.getHeight() > 1 && res.getWidth() > 1);
}

template<typename T>
size_t MipMap<T>::getOctaves() const{
    return _mipmap.size();
}
/**
 * @param size_des: provide desired size along dimension
        min_dim = min(ref_img.height, ref_img.width)
*/
template<typename T>
void MipMap<T>::get_interpolated_img(float ratio) const requires(std::is_same_v<T, PixelMap<png_byte>>){
    size_t height_des = ratio * _mipmap.at(0).getHeight();
    size_t width_des = ratio * _mipmap.at(0).getWidth();
    size_t shortest_dim = std::min(height_des, width_des);
    // recall smallest image at last position
    size_t offset = 1 + static_cast<size_t>(std::log2(static_cast<float>(shortest_dim)) / std::log2(2.0f));
    size_t idx_img = _mipmap.size()-offset;
    size_t idx_IMG = idx_img-1;
    printf("ref: (%ld, %ld) | target: (%ld, %ld) | REF: (%ld, %ld)\n", 
        _mipmap.at(idx_img).getHeight(), _mipmap.at(idx_img).getWidth(), 
        idx_img, idx_IMG,
        _mipmap.at(idx_IMG).getHeight(), _mipmap.at(idx_IMG).getWidth());
    //rescale::trilinear();
}

template<typename T>
void MipMap<T>::get_interpolated_img(float ratio) const requires(!std::is_same_v<T, PixelMap<png_byte>>){
    printf("wrong function pendejo\n");
    return;
}

template<typename T>
T& MipMap<T>::operator[](size_t idx){ return _mipmap[idx]; }

template<typename T>
const T& MipMap<T>::operator[](size_t idx) const{ return _mipmap[idx]; }

template<typename T>
void MipMap<T>::savefig(std::string_view filename) const requires(std::is_same_v<T, PixelMap<png_byte>>){

    size_t channels_tot = _mipmap.at(0).getChannels();
    size_t width_tot = 0;
    std::for_each(_mipmap.begin(), _mipmap.end(), [&](const auto& img){ width_tot += img.getWidth(); });
    size_t height_tot = _mipmap.at(0).getHeight() * channels_tot;
    
    // START: final image construction
    // height of the largest image in the mipmap
    size_t height_offset = _mipmap.at(0).getHeight();
    // need array of widths since offsets vary. Initialise first width offset
    std::vector<size_t> width_offset(_mipmap.size(), 0);
    // accumulate the sum of widths
    for(size_t n = 1; n < _mipmap.size(); ++n)
        width_offset[n] = width_offset[n-1] + _mipmap.at(n-1).getWidth();

    // final image
    PixelMap<png_byte> pixelmap_tot(height_tot, width_tot, channels_tot);
    
    // col. type [c=channel] (c=0:gray | c=1:gray+alpha | c=2:rgb | c=3:rgb+alpha)
    int color_type = png_utils.col_type[channels_tot-1];
    switch(color_type){
        case PNG_COLOR_TYPE_GRAY:
            printf("gray\n");
            // gray channel
            for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx)
                for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i, j+width_offset[img_idx], 0) = (_mipmap.at(img_idx))(i, j, 0);
            break;
        case PNG_COLOR_TYPE_GRAY_ALPHA:
            printf("gray + alpha\n");
            // gray channel
            for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx){
                for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i){
                    for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j){
                        pixelmap_tot(i, j+width_offset[img_idx], 0) = (_mipmap.at(img_idx))(i, j, 0);
                        pixelmap_tot(i, j+width_offset[img_idx], 1) = 0xff;
                    }
                }
            }
            // alpha channel
            for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx)
                for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i+height_offset, j+width_offset[img_idx], 1) = (_mipmap.at(img_idx))(i, j, 1);
            break;
        case PNG_COLOR_TYPE_RGB:
            printf("rgb\n");
            for(size_t c = 0; c < 3; ++c)
                for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx)
                    for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i)
                        for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j)
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], c) = (_mipmap.at(img_idx))(i, j, c);
            break;
        case PNG_COLOR_TYPE_RGB_ALPHA:
            printf("rgb + alpha\n");
            // rgb channels
            for(size_t c = 0; c < 3; ++c){
                for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx){
                    for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i){
                        for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j){
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], c) = (_mipmap.at(img_idx))(i, j, c);
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], 3) = 0xff;
                        }
                    }
                }
            }
            // alpha channel
            for(size_t img_idx = 0; img_idx < _mipmap.size(); ++img_idx)
                for(size_t i = 0; i < _mipmap.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < _mipmap.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i+3*height_offset, j+width_offset[img_idx], 3) = (_mipmap.at(img_idx))(i, j, 3);
            break;
        default:
            printf("unknown number of channels\n");
            exit(1);
    }
    // END: final image construction
    // Write image to png file
    write_png_file(filename, pixelmap_tot);
}

template<typename T>
void MipMap<T>::savefig(std::string_view filename) const requires(!std::is_same_v<T, PixelMap<png_byte>>){

    // START: BitMap -> PixelMap (conversion)
    std::deque<PixelMap<png_byte>> pixelmaps;
    pixelmaps.resize(_mipmap.size());
    auto it_bitmap = _mipmap.begin();
    auto it_pixelmap = pixelmaps.begin();
    for(; it_pixelmap != pixelmaps.end(); ++it_bitmap, ++it_pixelmap){
        *it_pixelmap = PixelMap<png_byte>(*it_bitmap); // think about making rvalue operator=
    }
    // END: BitMap -> PixelMap (conversion)

    size_t channels_tot = pixelmaps.at(0).getChannels();
    size_t width_tot = 0;
    std::for_each(_mipmap.begin(), _mipmap.end(), [&](const auto& img){ width_tot += img.getWidth(); });
    size_t height_tot = pixelmaps.at(0).getHeight() * channels_tot;
    
    // START: final image construction
    // height of the largest image in the mipmap
    size_t height_offset = pixelmaps.at(0).getHeight();
    // need array of widths since offsets vary. Initialise first width offset
    std::vector<size_t> width_offset(pixelmaps.size(), 0);
    // accumulate the sum of widths
    for(size_t n = 1; n < pixelmaps.size(); ++n)
        width_offset[n] = width_offset[n-1] + pixelmaps.at(n-1).getWidth();

    // final image
    PixelMap<png_byte> pixelmap_tot(height_tot, width_tot, channels_tot);
    
    // col. type [c=channel] (c=0:gray | c=1:gray+alpha | c=2:rgb | c=3:rgb+alpha)
    int color_type = png_utils.col_type[channels_tot-1];
    switch(color_type){
        case PNG_COLOR_TYPE_GRAY:
            printf("gray\n");
            // gray channel
            for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx)
                for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i, j+width_offset[img_idx], 0) = (pixelmaps.at(img_idx))(i, j, 0);
            break;
        case PNG_COLOR_TYPE_GRAY_ALPHA:
            printf("gray + alpha\n");
            // gray channel
            for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx){
                for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i){
                    for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j){
                        pixelmap_tot(i, j+width_offset[img_idx], 0) = (pixelmaps.at(img_idx))(i, j, 0);
                        pixelmap_tot(i, j+width_offset[img_idx], 1) = 0xff;
                    }
                }
            }
            // alpha channel
            for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx)
                for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i+height_offset, j+width_offset[img_idx], 1) = (pixelmaps.at(img_idx))(i, j, 1);
            break;
        case PNG_COLOR_TYPE_RGB:
            printf("rgb\n");
            for(size_t c = 0; c < 3; ++c)
                for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx)
                    for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i)
                        for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j)
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], c) = (pixelmaps.at(img_idx))(i, j, c);
            break;
        case PNG_COLOR_TYPE_RGB_ALPHA:
            printf("rgb + alpha\n");
            // rgb channels
            for(size_t c = 0; c < 3; ++c){
                for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx){
                    for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i){
                        for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j){
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], c) = (pixelmaps.at(img_idx))(i, j, c);
                            pixelmap_tot(i+c*height_offset, j+width_offset[img_idx], 3) = 0xff;
                        }
                    }
                }
            }
            // alpha channel
            for(size_t img_idx = 0; img_idx < pixelmaps.size(); ++img_idx)
                for(size_t i = 0; i < pixelmaps.at(img_idx).getHeight(); ++i)
                    for(size_t j = 0; j < pixelmaps.at(img_idx).getWidth(); ++j)
                        pixelmap_tot(i+3*height_offset, j+width_offset[img_idx], 3) = (pixelmaps.at(img_idx))(i, j, 3);
            break;
        default:
            printf("unknown number of channels\n");
            exit(1);
    }
    // END: final image construction
    // Write image to png file
    write_png_file(filename, pixelmap_tot);
}