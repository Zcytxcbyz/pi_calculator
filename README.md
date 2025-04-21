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
     sudo apt-get install libgmp-dev cmake
     ```

   - On macOS:
     ```bash
     brew install gmp cmake
     ```

   - On Windows (using MSYS2):
     ```bash
     pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-cmake
     ```

2. **Clone the Repository**
   ```bash
   git clone <repository-url>
   cd pi_calculator
   ```

3. **Generate Build Files**
    ```bash
    cmake -S . -B build
    ```

4. **Build the Project**
    ```bash
    cmake --build build
    ```

5. **Run the Program**
    ```bash
    ./build/pi_calculator
    ```

## Usage

Run the program with the following options:

```bash
./pi_calculator [options]
```

### Options

- `-d <digits>`: Specify the number of digits to calculate (default: 1000).

- `-o <filename>`: Specify the output file name (default: pi.txt).

- `-t <threads>`: Specify the number of threads to use (default: number of CPU cores).

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

## Notes

- The higher the number of digits, the more memory and computation time required.

- The output file contains π formatted with 100 digits per line for readability.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
