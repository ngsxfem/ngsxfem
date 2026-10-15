from pathlib import Path

def get_cmake_dir():
    return str(Path(__file__).resolve().parent / "cmake")

def get_include_dir():
    return str(Path(__file__).resolve().parent / "include")

def get_library_dir():
    return str(Path(__file__).resolve().parent)

if __name__ == "__main__":
    print(get_cmake_dir())
