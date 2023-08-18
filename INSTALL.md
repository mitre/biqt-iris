# Linux (Ubuntu 22.04)

## Prerequisites

This guide targets Ubuntu 22.04

This provider relies on OpenCV and other core development tools including gcc. These dependencies can be installed from 
the default apt repository using the following commands:

TODO: Update

## Building and Installing

TODO: Update

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
