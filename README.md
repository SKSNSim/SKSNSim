# SKSNSim: vector generator of Supernova burst and Supernovae Diffused Neutrino Background
## Usage
There are two (four at this moment) binaries, "main\_snburst(\_new)" and "main\_dsnb(\_new)" in this package. They are for Supernova burst (SN burst) and for Supernovae Diffused Neutrino Backgroud (DSNB) individually. The suffix of "\_new" means newly implemented version. We recoomend to use the new version.

### Compile
It is based on GNUMake. Just executing ``make`` command. However, it requires SK offline software packages.
```SHELL
$ make clean # for clean-up
$ make # for compiling
```
Compiled binaries are places into ``bin`` directory.

### main\_snburst(\_new)
For previous version, ``main\_snburst``, the basic usage is ``main\_snburst {model} {nuosc} {distance} {generatormode} {outputdirectory} {randomseed}``. List of arguments are below:
* 1st: model name (e.g. nakazato/intp2001.data)
* 2nd: neutrino-oscillation (0:no, 1:normal, 2:inverted)
* 3rd: Distance (normalized to 10kpc)
* 4th: Generate event or not (0: No just calculate expected, 1: Yes for detector simulation)
* 5th: Output directory, default: ./data/(SN model)
* 6th: Random seed (not supported yet)

For new version, you can see detail usage via executing ``main\_snburst\_new --help``. Basically, a format of the previous version is available. But we recommend to use new format which can be specified via options starting ``--`` (e.g. ``--seed``). All options have default value, and you do not have to specify all of them.


### main\_dsnb(\_new)
Detail can be dumped by executing ``main\_dsnb(\_new) --help``.
Basically, you can run by just executing ``main\_dsnb(\_new)`` with wanted options you would like to change.
