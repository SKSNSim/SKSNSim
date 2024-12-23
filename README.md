# SKSNSim: vector generator of Supernova burst and Diffuse Supernova Neutrino Background
## Usage
There are two binaries, ``main_snburst`` and ``main_dsnb`` in this package. They are for Supernova burst (SN burst) and for Diffuse Supernova Neutrino Backgroud (DSNB) individually.

### Compile
It is based on GNUMake. Just executing ``make`` command. However, it requires SK offline software packages.
```SHELL
$ make clean # for clean-up
$ make # for compiling
```
Compiled binaries are places into ``bin`` directory.

### main\_snburst

You can see detail usage via executing ``main_snburst --help``. All options have default value, and you do not have to specify all of them.


### main\_dsnb
Detail can be dumped by executing ``main_dsnb --help``.
Basically, you can run by just executing ``main_dsnb`` with wanted options you would like to change.

## More detail:

The guide for users and developpers are available on https://github.com/SKSNSim/SKSNSim/releases/download/v1.2.0/guide_sksnsim.pdf , which is output of doc/guide_sksnsim.texi (Texinfo file). If you want to see the PDF version, please do ``make doc``.

You can also see the SKSNSim paper from https://iopscience.iop.org/article/10.3847/1538-4357/ad344e.
