# aaHMM

`alpha-HMM` is a tool that RNA secondary structure prediciton.

## Requirements

- C++ with C++17 support
- CMake version >= 3.26.5
- Python 3.8 or newer (with `subprocess` module)


## Building the C++ Tool

1. From the project root directory, navigate to the `aaHMM` directory:
   ```bash
   cd aaHMM
2. Create a build directory and navigate into it:

    ```bash
    mkdir build
    cd build
3. Run CMake and then compile:

    ```bash
    cmake ..
    make
4. The executable `aaHMM` should now be present in the `build` directory.

## Using the Python Wrapper

We provide a Python script that interfaces with the `alphaHMM` executable to process sequences and obtain predictions.

### Example

To predict the structure for a given sequence, use the following code:

    
    import subprocess
    
    executable_path = "aaHMM /cmake-build-debug/aaHMM"
    
    def predict_from_cpp(sequence):
        result = subprocess.run([executable_path, sequence], capture_output=True, text=True)
        return result.stdout.strip()
    
    sequence = "UAUAUUUUAGUGUAUGAUGCACAAAAGUUUUUGAAACUUUUAGAAAUAGUUUAAUUCUAUUAAAUAUAACCA"
    prediction = predict_from_cpp(sequence)
    
    print("Prediction:")
    print(prediction)
    

Ensure the `alphaHMM` executable is located as per the `executable_path` or adjust the path accordingly.


## Contact & Support

For any issues or queries, contact [Sixiang] at [sz94706@uga.edu].
