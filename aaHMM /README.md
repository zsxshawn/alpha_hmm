# aaHMM

`aaHMM` is a tool that [BRIEF DESCRIPTION OF WHAT YOUR TOOL DOES].

## Requirements

- C++ with C++17 support
- CMake version >= 3.26.5

## Building the C++ Tool

1. Navigate to the project root directory.
2. Create a build directory and navigate into it:

\```bash
mkdir build
cd build
\```

3. Run CMake and then compile:

\```bash
cmake ..
make
\```

4. The executable `aaHMM` should now be present in the `build` directory.

## Using the Python Wrapper

We provide a Python script that interfaces with the `aaHMM` executable to process sequences and obtain predictions.

### Example

To predict the structure for a given sequence, use the following code:

\```python
import subprocess

executable_path = "./build/aaHMM"

def predict_from_cpp(sequence):
result = subprocess.run([executable_path, sequence], capture_output=True, text=True)
return result.stdout.strip()

sequence = "UAUAUUUUAGUGUAUGAUGCACAAAAGUUUUUGAAACUUUUAGAAAUAGUUUAAUUCUAUUAAAUAUAACCA"
prediction = predict_from_cpp(sequence)

print("Prediction:")
print(prediction)
\```

Ensure the `aaHMM` executable is located as per the `executable_path` or adjust the path accordingly.

## License

[Information about the project's license.]

## Contact & Support

For any issues or queries, contact [Your Name] at [Your Email Address].
