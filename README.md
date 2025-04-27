# Pi Calculator

`pi_calculator` is a high-precision tool for calculating the value of π (pi) using the Chudnovsky algorithm. It leverages the GNU MP (GMP) library for arbitrary-precision arithmetic and supports multithreaded computation for improved performance.

## Features

- Calculates π to a specified number of digits using the Chudnovsky algorithm.
- Supports multithreaded computation with OpenMP.
- Outputs the result to a file in a formatted manner.
- Highly efficient and scalable for large computations.

## Dependencies

- [GMP (GNU Multiple Precision Arithmetic Library)](https://gmplib.org/)
- [OpenMP](https://www.openmp.org/) (for multithreading support)
- A C compiler supporting the C11 standard.

## Build Instructions

### Using CMake

1. **Install Dependencies**

   - On Linux:
     ```bash
     sudo apt-get install libgmp-dev cmake build-essential
     ```

   - On macOS:
     ```bash
     brew install gmp cmake
     ```

   - On Windows (using MSYS2):
     ```bash
     pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-cmake
     ```

   - On Windows (using Visual Studio):
     1. Install [Visual Studio](https://visualstudio.microsoft.com/) with the "Desktop development with C++" workload.
     2. Install [CMake](https://cmake.org/download/).
     3. Download and build the GMP library manually or use a precompiled version.

2. **Clone the Repository**
   ```bash
   git clone https://github.com/Zcytxcbyz/pi_calculator.git
   cd pi_calculator
   ```

3. **Generate Build Files**
   - On Linux/macOS:
     ```bash
     cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON
     ```

   - On Windows (using MSYS2):
     ```bash
     cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON
     ```

   - On Windows (using Visual Studio):
     ```bash
     cmake -S . -B build -G "Visual Studio 16 2019" -A x64 -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON
     ```

4. **Build the Project**
   - On Linux/macOS:
     ```bash
     cmake --build build
     ```

   - On Windows (using MSYS2):
     ```bash
     cmake --build build
     ```

   - On Windows (using Visual Studio):
     ```bash
     cmake --build build --config Release
     ```

5. **Run the Program**
   - On Linux/macOS:
     ```bash
     ./build/pi_calculator
     ```

   - On Windows:
     ```bash
     .\build\Release\pi_calculator.exe
     ```

---

### Notes

- **Static Linking**:
  - By default, the project is configured to build a static executable (`BUILD_STATIC=ON`).
  - On Windows with MSVC, the `/MT` flag is used for static runtime linking.
  - On Linux/macOS, `-static`, `-static-libgcc`, and `-static-libstdc++` are used for static linking.

- **Debug Build**:
  - To build in debug mode, replace `-DCMAKE_BUILD_TYPE=Release` with `-DCMAKE_BUILD_TYPE=Debug`.

- **Dependencies**:
  - Ensure that the GMP library is correctly installed and accessible by the compiler.
  - OpenMP support is required for multithreaded computation.

## Usage

Run the program with the following options:

```bash
./pi_calculator [options]
```

### Options

- `-d <digits>`: Specify the number of digits to calculate (default: 1000).

- `-o <filename>`: Specify the output file name (default: pi.txt).

- `-t <threads>`: Specify the number of threads to use (default: number of CPU cores).

- `-f`: Format output (default: unformatted);

- `-c`: Disable output file

- `-b <size>`: Set buffer size in bytes (default: 65536)

- `-s <schedule>`: Set OpenMP schedule type (static, dynamic, guided) and chunk size (default: guided)

- `-h`: Display the help message.

### Examples

1. Calculate π to 5000 digits and save the result to `output.txt`:
    ```bash
    ./pi_calculator -d 5000 -o output.txt
    ```

2. Use 8 threads to calculate π to 10000 digits:
    ```bash
    ./pi_calculator -d 10000 -t 8
    ```

3. Calculate π to 1000 digits without saving to a file:
    ```bash
    ./pi_calculator -d 1000 -c
    ```

## Performance Notes

- The higher the number of digits, the more memory and computation time required.

- Multithreading can significantly speed up calculations, especially on systems with multiple CPU cores.

- The caching mechanism (`ENABLE_CACHE`) can optimize repeated calculations for large values of `k`.

## Build Options

The project supports several build options that can be configured using CMake:

- `ENABLE_SIMD`: Enable SIMD instructions (e.g., AVX2/SSE4) for performance optimization (default: ON).

- `ENABLE_LTO`: Enable Link Time Optimization (LTO) for smaller and faster binaries (default: ON).

- `BUILD_STATIC`: Build as a statically linked executable (default: ON).

- `ENABLE_CACHE`: Enable caching for large calculations (default: ON).

- `OUTPUT_FORMAT`: Enable formatted output (e.g., 100 digits per line) (default: OFF).

- `MULTI_OUTPUT`: Enable multi-output support (default: ON).

To enable or disable these options, pass `-D<option>=ON/OFF` to the `cmake` command. For example:

```bash
cmake -S . -B build -DENABLE_SIMD=OFF -DBUILD_STATIC=OFF
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.
