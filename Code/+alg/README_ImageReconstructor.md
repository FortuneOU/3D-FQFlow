## 🖼️ Image Reconstruction & Processing Module

Reconstructing three-dimensional ultrafast Power Doppler (uPDI) sequences from massive raw Radio-Frequency (RF) data is notoriously time-consuming. Traditional Delay-and-Sum (DAS) beamforming recalculates time-of-flight delays for every pixel/voxel across hundreds of frames, which can take days for high-resolution 3D volumes. Our framework overcomes this bottleneck by implementing a highly optimized, matrix-based DAS acceleration strategy combined with advanced post-processing tools.

### ✨ Key Features & Technical Optimizations

1. **Precomputed Sparse DAS Matrix Acceleration:** Based on the linear superposition principle, the beamforming process for any spatial location is independent. Instead of loop-based recalculations, the system precomputes a massive, sparse mapping matrix for each transmit angle. 
2. **Ultra-Fast Frame Processing:** Once the delay-and-sum matrix is stored in memory, newly simulated RF or I/Q frames are beamformed using a single, highly optimized sparse matrix multiplication. This reduces the reconstruction time of a 30-frame 3D sequence from over 16 hours down to under an hour.
3. **Comprehensive Post-Processing & Evaluation:** The module natively supports Singular Value Decomposition (SVD) clutter filtering to separate blood flow from tissue motion, alongside B-mode logarithmic compression. It also includes built-in quantitative evaluation tools providing MSE, PSNR, and SSIM metrics against the hemodynamic ground truth.

### 💻 Core Code Parsing (`FQdasmtx3` Function)

The `FQdasmtx3` function is the architectural core of the reconstruction module. It is responsible for generating the highly efficient Delay-and-Sum mapping matrix (`M`) for 3D imaging using matrix arrays. Its workflow is structured as follows:

1. **Spatial & Parameter Initialization:** The function ingests the target 3D spatial grid coordinates `(X, Y, Z)`, transducer element properties, sampling frequencies, and transmit delay laws. It handles both raw RF signals and complex I/Q data (automatically applying necessary phase rotations for the latter).
2. **Time-of-Flight & Sub-aperture Filtering:** It calculates the precise transmit and receive distances (`dTX` and `dRX`) for every element-to-voxel pair. It then converts these acoustic travel times into fast-time indices. Crucially, it applies dynamic receive aperture limits (using the `fnumber` parameter) to exclude signals outside the valid directivity angles, enhancing image resolution and reducing matrix density.
3. **Multi-Method Sub-sample Interpolation:** To ensure phase accuracy without oversampling the raw data, the function maps fractional time indices using user-defined interpolation kernels. It supports methods ranging from basic 'nearest' and 'linear' (fastest) to high-fidelity 'quadratic', '5points', and windowed sinc functions like 'lanczos3' and 'lanczos5'.
4. **Sparse Matrix Assembly & Execution:** The computed interpolation weights and indices are assembled into an enormous but highly memory-efficient sparse matrix `M`. In the main reconstruction loop, the entire 3D volume for a given frame is instantly formed by executing the algebraic operation `bfSIG = M * SIG(:)`, dramatically shifting the computational load from iterative logic to optimized linear algebra.
