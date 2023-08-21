# Linux (Ubuntu 22.04)

## Prerequisites

This guide targets Ubuntu 22.04

This provider relies on OpenCV and other core development tools including gcc. These dependencies can be installed using the following commands:

```bash
sudo apt update
sudo apt install -y cmake build-essential libopencv-dev
```

## Building and Installing

Clone the repository.

```bash
git clone https://github.com/mitre/biqt-iris.git
cd biqt-iris
```

Finally, build and install the provider. Set BIQT_HOME via environment variable or using -DBIQT_HOME in cmake command below.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
```

## Running the Provider

After installation, you can call the provider using the reference BIQT CLI as follows:

```
# Runs all iris providers (including this one).
$> biqt -m "iris" <image path>

# Runs only the biqt-iris provider.
$> biqt -p "BIQTIris" <image path>
```

# Windows

Windows is no longer a supported platform.
