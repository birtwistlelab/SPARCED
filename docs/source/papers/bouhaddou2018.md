# Bouhaddou 2018 model

Files to exactly replicate Bouhaddou2018 model are located in SPARCED Brep
[TODO: explain how to get back to prior version].

In this format, mRNAs and protein species are tracked in separate variables.
The number of species and ratelaws are also the same Bouhaddou2018 model.
SPARCED is the new, updated and cleaned up version.

## Model Testing and Performance

The Docker image of the model and simulations are tested with multiple machines
(see bellow). The simulation times are:

```
- Model file creation time: ~ 2 min
- Model compilation time: 15-25 min (depends on model size)
- Model simulation time: ~ 1 min for 24-hour simulation
```

1. Ubuntu-Desktops:
  - Ubuntu 18.04, Intel Core i7 3930 CPU @ 3.20 Ghz, 32 GB, DDR3, Nvidia GTX
  690 GPU
2. Windows-Desktops:
  - Windows 10 Education, Intel Core i5-3470 CPU @ 3.20 Ghz, 8.00 GB RAM,
  Nvidia GTX 650 GPU, 64-bit operating system
  - Windows 10 Education, Intel Core i7-9700 CPU @5.00 GHz, 32.00 GB RAM,
  Nvidia RTW 2070 GPU, 64-bit operating system
3. Windows-Laptops:
  - Windows 10 Pro, Intel Core i7-8550U CPU @ 2.00 GHz, 16.00 GB RAM, Intel UHD
  620 GPU, 64-bit operating system
  - Windows 10 Education, Intel Core i7-10705H CPU @ 2.60 GHz, 16.00 GB RAM,
  Nvidia GTX 1660ti GPU, 64-bit operating system

