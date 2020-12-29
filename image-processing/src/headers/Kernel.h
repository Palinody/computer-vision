#pragma once

//#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>

// non random distributions
enum Distr : unsigned char{ GAUSS };

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
//////////////////////////////////////////////////////
//                  kernel mask                     //
//////////////////////////////////////////////////////

template<bool Buf = true>
class ker_mask{
public:
    explicit
    ker_mask(size_t k_size, size_t k_idx, size_t size) requires(!Buf);
    ker_mask(size_t k_size, size_t k_idx, size_t size) requires(Buf);

    inline void update(size_t img_idx);
    inline size_t& operator()(size_t k) requires(!Buf);
    inline const size_t& operator()(size_t k) const requires(!Buf);
	inline size_t& operator()(size_t img_idx, size_t k) requires(Buf);
    inline const size_t& operator()(size_t img_idx, size_t k) const requires(Buf);

private:
    // kernel size, reference ("middle") index in kernel
    size_t _k_size, _k_idx;
    // convolved image size
    size_t _size;
    // number of masks in case of buffered class (else 1)
    size_t _masks;
    // reference indices lookup table
    std::unique_ptr<size_t[]> _ker_indices;
    // runtime indices lookup table updated 
    // w.r.t. curr. image idx or buffered 
    // and accessed w.r.t. curr img idx
    std::unique_ptr<size_t[]> _ker_indices_update;
};

template<bool Buf>
ker_mask<Buf>::ker_mask(size_t k_size, size_t k_idx, size_t size) requires(!Buf) :
        _k_size{ k_size }, _k_idx{ k_idx }, _size{ size }, _masks{ 1 },
        _ker_indices{ std::make_unique<size_t[]>(k_size) },
        _ker_indices_update{ std::make_unique<size_t[]>(k_size) }{
    
    for(size_t k = 0; k < _k_idx; ++k) _ker_indices[k] = _k_idx - k;
    for(size_t k = _k_idx; k < _k_size; ++k) _ker_indices[k] = k - _k_idx;
}

template<bool Buf>
ker_mask<Buf>::ker_mask(size_t k_size, size_t k_idx, size_t size) requires(Buf) : 
        _k_size{ k_size }, _k_idx{ k_idx }, _size{ size }, _masks{ std::min(k_size, size) },
        _ker_indices{ std::make_unique<size_t[]>(k_size) },
        _ker_indices_update{ std::make_unique<size_t[]>(k_size*_masks) }{
    
    for(size_t k = 0; k < _k_idx; ++k) _ker_indices[k] = _k_idx - k;
    for(size_t k = _k_idx; k < _k_size; ++k) _ker_indices[k] = k - _k_idx;
    
    for(size_t m = 0; m < _masks; ++m){
        for(size_t k = 0; k < _k_idx; ++k)
			_ker_indices_update[k+m*_k_size] = std::min(_ker_indices[k], m);
        for(size_t k = _k_idx; k < _k_size; ++k)
			_ker_indices_update[k+m*_k_size] = std::min(_ker_indices[k], _masks-1-m);
    }
}

template<bool Buf>
void ker_mask<Buf>::update(size_t img_idx){
    
    for(size_t k = 0; k < _k_idx; ++k)
        _ker_indices_update[k] = std::min(_ker_indices[k], img_idx);
    for(size_t k = _k_idx; k < _k_size; ++k)
        _ker_indices_update[k] = std::min(_ker_indices[k], _size-1-img_idx);
}

template<bool Buf>
size_t& ker_mask<Buf>::operator()(size_t k) requires(!Buf){
    
    return _ker_indices_update[k];
}

template<bool Buf>
const size_t& ker_mask<Buf>::operator()(size_t k) const requires(!Buf){
    
    return _ker_indices_update[k];
}

template<bool Buf>
size_t& ker_mask<Buf>::operator()(size_t img_idx, size_t k) requires(Buf){

	size_t m = (img_idx <= _k_idx) ? 
        img_idx : ((_size-1-img_idx < _k_size-_k_idx) ? 
            _masks - (_size-img_idx) : _k_idx);
    return _ker_indices_update[k+m*_k_size];
}

template<bool Buf>
const size_t& ker_mask<Buf>::operator()(size_t img_idx, size_t k) const requires(Buf){

    size_t m = (img_idx <= _k_idx) ? 
        img_idx : ((_size-1-img_idx < _k_size-_k_idx) ? 
            _masks - (_size-img_idx) : _k_idx);
    return _ker_indices_update[k+m*_k_size];
}