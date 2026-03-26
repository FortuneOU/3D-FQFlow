## 📡 Ultrasound Simulation Module

One of the primary bottlenecks in 3D ultrafast Power Doppler Imaging (uPDI) simulation is the enormous computational cost and memory consumption required to simulate hundreds of frames for complex 3D structures. To overcome this limitation, our framework introduces a highly optimized, GPU-accelerated frequency-domain ultrasound simulator based on linear acoustic theory and the weak scattering assumption.

### ✨ Key Features & Technical Optimizations

1. **GPU-Accelerated Acoustic Solver:** The core simulator is built upon an optimized PFIELD approach. The most computationally intensive operations—specifically, the matrix multiplications between propagation accumulation terms and incremental frequency-stepping terms—are entirely offloaded to the GPU. This reduces the time required to simulate millions of scatterers from hours to just a few seconds.
2. **Dynamic Memory Partitioning:** Traditional simulators often suffer from Out-Of-Memory (OOM) failures when processing large 3D volumes. To prevent this, our simulator dynamically detects available GPU memory and partitions the massive scatterer arrays into optimal, bite-sized chunks for distributed parallel computation.
3. **Multi-Plane Transmission (MPT):** The simulator natively supports multi-angle plane-wave transmissions required for ultrafast imaging, efficiently processing complex transmit delay laws and reconstructing the corresponding radio-frequency (RF) fields.

### 💻 Core Code Parsing (`FQsimus3` Function)

The `FQsimus3` function serves as the central engine for generating 3D RF signals. It simulates the echoes generated when a planar 2D array transmits specified wavefronts through a scattering medium. Its execution pipeline is structured as follows:

1. **Input Parsing & Setup:** The function receives scatterer spatial coordinates `(X, Y, Z)`, reflection coefficients `(RC)`, transmit delays, and transducer parameters. It automatically calculates the maximum acoustic propagation depth and determines the optimal frequency sampling step (`df`) to strictly avoid aliasing in the time domain.
2. **GPU Chunking & Distributed Computation:** 
   - The script first queries the hardware capabilities using `gpuDevice.AvailableMemory`.
   - It then computes the required memory footprint per scatterer (`MemoryUnit`).
   - Based on the total scatterer count (`Nx`), the data is automatically divided into `NW` distributed chunks (`numPointsCal`).
   - A loop sequentially feeds these data chunks into the underlying solver (`alg.FQpfield3`), which accumulates the RF spectrum (`RFspectrum`) piece by piece. This guarantees stable execution even on consumer-grade GPUs with limited VRAM.
3. **Time-Domain Reconstruction & Denoising:** Once the full frequency-domain spectrum is calculated, the function performs an Inverse Fast Fourier Transform (`ifft`) with a symmetric flag to accurately reconstruct the time-domain RF data matrix. Finally, a nonlinear hyperbolic tangent (`tanh`) thresholding function is applied to gracefully suppress near-zero numerical noise.
