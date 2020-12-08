> C++ repository used for educational purposes in computer vision. Contains a variety of algorithms that will be used for a SLAM project using [SIFT](https://en.wikipedia.org/wiki/Scale-invariant_feature_transform) algorithm. The ultimate goal of the project is to launch it on embedded systems (smartphones, robots)

### TODO
fix: clang++

fix: fft works with images with limited size on my machine: approx. 900x900x4 (no memory optim.)

add: gray and gray + alpha bitmaps (rest of lib already took them into consideration, no mod. needed)

do: all the rest...


### Compilation
```sh
$ g++-10 -std=c++20 -Wall -O3 -march=native -o main main.cpp
$ ./a.out
```

### Requirements
- c++-20
- g++-10 (recom. bec. only compiler tested and we need c++-20)
- [libpng](https://github.com/glennrp/libpng)

### Custom containers
- PixelMap<png_byte>: color channels read/wrote via indexing
- BitMapRGBA: color channels read/wrote via bit shifting (channels fixed to 4, type fixed to uint32_t)
- BitMapRGB: color channels read/wrote via bit shifting (channels fixed to 3, type fixed to uint32_t)
- MipMap: Pyramid of previous image containers

### Streams
- Reads/writes from/to png files with libpng library and custom containers

### Image processing algorithms
#### Rescale (any scale)
- Nearest neighbor
- Bilinear
- Trilinear
- Bicubic

#### Image processing
- Fast Fourier Transform (FFT) with gaussian kernel
- Difference of Gaussians (DoG) with FFT

### Some results
#### MipMap
![mip_map](https://user-images.githubusercontent.com/32341154/101549555-3028b900-39ae-11eb-95d8-7f0bf5b37c22.png)

#### Difference of Gaussians (DoG)
##### Original
![000369_left](https://user-images.githubusercontent.com/32341154/101549849-ac230100-39ae-11eb-9ae6-4b51f5591d78.png)
##### Processed
![diff_of_gaussians_rgb](https://user-images.githubusercontent.com/32341154/101549572-374fc700-39ae-11eb-8085-4a840dfa41cf.png)

#### Noisy Gaussian distribution
![gauss_ball](https://user-images.githubusercontent.com/32341154/101549549-2e5ef580-39ae-11eb-9f09-488ad09dc6ff.png)
