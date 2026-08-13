# Project Context: Enzo OpenMP GPU Offloading & Asynchronous RT

## Overview
This project extends Enzo's GPU compute capabilities using OpenMP target offloading and optimizes the Moray ray tracer. Key objectives include porting C++ hydro solvers to GPUs and restructuring the ray tracer for many-core efficiency and communication overlap.

## Key Source Files
- `src/`: Core solvers.
- `src/Grid_TransportPhotonPackages.C`: Primary ray propagation logic (SoA Target).
- `src/EvolvePhotons.C`: High-level RT coordination and MPI logic.
- `src/Grid_SolvePPM_DE.C`: Reference for stencil parallelization.

## Architecture Guidelines
- **Structure of Arrays (SoA):** Convert the photon package linked list into a flat SoA to enable coalesced memory access and SIMT parallelization [1, 2].
- **Communication-Computation Overlap:** Utilize an asynchronous, six-step MPI handshake to hide latency when rays cross AMR grid boundaries [3, 4].
- **Data Minimization:** Keep fields (Density, Energy, etc.) on-device during the RT step. Accumulate photo-ionization ($k_{ph}$) and heating ($\Gamma_{ph}$) rates in device-side buffers [5, 6].

## Development Rules
1. **SoA Layout:** Maintain the 11 essential fields: Flux, Type, Energy, Column Density, Pixel Num, Level (Single/Int) and TimeInterval, EmissionTime, CurrentTime, Radius, SourcePosition (Double/Precision-dependent) [7].
2. **Pre-allocation:** Account for ray splitting (1 ray to 4 child rays) by pre-allocating device buffers with sufficient capacity to avoid mid-step reallocations [8, 9].
3. **MPI Handshake:** 
   - Ranks must first exchange `Nmesg` (count of expected data packets) via non-blocking `MPI_Isend`/`MPI_Irecv` [4, 10].
   - Aggressively drain the message stack before checking data status [10].
4. **Asynchronous Loop:** Use `MPI_Testsome` while local rays are available for tracing; fall back to `MPI_Waitsome` only when local work is exhausted [11].
5. **Termination:** Ensure global termination checks include both empty local queues and zero outstanding MPI requests to prevent race conditions [12].

## Primary Targets
1. **PPM/ZEUS Offloading:** Port 1D sweep loops and staggered mesh updates to OpenMP target regions.
2. **Ray Tracer SoA:** Replace the `PhotonPackage` linked list in `Grid_TransportPhotonPackages.C` with the SoA model for GPU offloading.
3. **Non-blocking MPI:** Implement the 6-step asynchronous communication protocol in `EvolvePhotons.C` to eliminate blocking boundary synchronization.