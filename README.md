## Linux Build Steps
This assumes available Netgen/NGSolve and Boost-libraries

```
git clone https://gitlab.asc.tuwien.ac.at/jschoeberl/xfem.git
cd amg
mkdir build
cd build
cmake ../ && make && make install
```

To build without tests in the last step use: `cmake -DINSTALL_DIR=../xfem-inst ../xfem-src`

## Examples
To run the python examples be sure to follow the build steps above.
Then navigate into the `py_tutorials` or `stokes/py_demos` or `tracefem/py_demos` and run
`netgen example.py`
where `example.py` stands for any of the available python files.

## Testing
Tests are enabled by default.
To run the test navigate to the build directory and run `make test`
or `ctest`.
If you need to see specific tests failing use `ctest -V`.
To run individual tests use `ctest -R <regex>`. E.g. `ctest -R h1` to only run h1 integration
tests.
