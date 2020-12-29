#pragma once

//#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>

// non random distributions
enum Distr{ GAUSS };

template<typename T>
class Kernel1D{
public:
    Kernel1D(size_t length, size_t channels);
    Kernel1D(size_t length, size_t channels, T sd);

    inline size_t getLength() const{ return _length; }
    inline size_t getChannels() const{ return _channels; }

    inline T& operator()(size_t idx, size_t c){ return _kernel[c + idx * _channels]; }
    inline const T& operator()(size_t idx, size_t c) const{ return _kernel[c + idx * _channels]; }

    void print() const;

private:
    size_t _length, _channels;
    std::unique_ptr<T[]> _kernel;
};

template<typename T>
Kernel1D<T>::Kernel1D(size_t length, size_t channels) : 
        _length{ length }, _channels{channels},
        _kernel{ std::make_unique<T[]>(length*channels) }{}

// mean of gauss. distr. always at the center of the ker.
template <typename T>
Kernel1D<T>::Kernel1D(size_t length, size_t channels, T sd) : 
    _length{length}, _channels{channels},
    _kernel{std::make_unique<T[]>(length * channels)}{

    // 2*sigma²
    T tss = 2 * sd * sd;
    // 1 / sqrt(2*pi*sd*sd)
    T norm = 1 / std::sqrt(static_cast<T>(M_PI) * tss);
    tss = 1 / tss;
    // if pixels even then the mean lies "between" 2 pixels
    T offset = (_length % 2 == 0) ? 0.5 : 0;
    // origin shifted by _pixels / 2.0f
    T shift = _length / 2 - offset;

    T sum = 0;
    for (size_t idx = 0; idx < _length; ++idx){
        T val = norm * std::exp(-(idx - shift) * (idx - shift) * tss);
        for (size_t c = 0; c < channels; ++c){
            _kernel[c + idx * _channels] = val;
        }
        sum += val;
    }
    // normalizing so that sum(ker) = 1
    for (size_t n = 0; n < _length * _channels; ++n)
        _kernel[n] /= sum;
}

template<typename T>
void Kernel1D<T>::print() const{

    for(size_t idx = 0; idx < _length; ++idx){
        printf("{");
        for(size_t c = 0; c < _channels; ++c){
            printf("%.8f, ", static_cast<float>(_kernel[idx * _channels]));
        }
        printf("},\n");
    }
    printf("\n");
}