import subprocess

executable_path = "aaHMM /cmake-build-debug/aaHMM"

def predict_from_cpp(sequence):
    result = subprocess.run([executable_path, sequence], capture_output=True, text=True)
    return result.stdout.strip()

sequence = "UAUAUUUUAGUGUAUGAUGCACAAAAGUUUUUGAAACUUUUAGAAAUAGUUUAAUUCUAUUAAAUAUAACCA"
prediction = predict_from_cpp(sequence)

print("Prediction:")
print(prediction)
