# Thesis2022
Thesis work for bachelor degree - Modeling Lunar Gravitational Perturbations in Low Orbit Integration

This repo consists of a pdf file providing the final thesis paper presented at the graduation commission - Be aware that some translating errors may be present. Also, we are required to not go beyond the 7 page length by the university itself, so it was not a lack of effort. I highly suggest to start from here


The Serial and Parallel folder contains the implementation of the Pines algorithm on the CPU and GPU (using CUDA) respectively, while the Disk and ModifiedDisk folder contains the main classes and methods that differentiates them from the basic Pines implementation.
Orbit propagation with the Disk and ModifiedDisk algorithms were implemented by modifying accordingly the source code of GMAT, but i lost the modified files to do a diff with, so you just have to trust me on that or you can modify your GMAT source as well - hint: interested GMAT files are at /src/base/forcemodel/harmonic/ of the main repo.
