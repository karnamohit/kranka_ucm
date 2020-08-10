# kranka_ucm
Codes written at UC Merced by and for the Isborn research group. Most are compatible with `python3.x`
- `scripts/`: these mostly interface with the `Gaussian` program. There are some standalone codes for running TDCI calculations as well.
- `ML/`: these codes are written using `tensorflow` and `pytorch` libraries (to simulate deep, feed-forward static neural networks for time-series analysis). Inspired and modified by [`@hbhat4000`](https://github.com/hbhat4000)
- `examples/`: this directory contains example computations run using the scripts and codes stored in the above folder(s). `examples/TDCI` refers to an example of a TDCI on the LiH molecule using the STO-3G basis set. This starts with using the configuration interaction scripts, `calc_tdm_cis_latest.ipynb` and `calc_tdm_casscf_latest.ipynb`, to extract data from Gaussian DV calculation `.log` files and storing the data in `tdm_tdci.txt` and CI density matrices in a series of `.npz` files, and then running the time-dependent  configuration interaction script, `TDCI_1pulse_latest_compact.ipynb`, which reads `tdm_tdci.txt` and using the CI density matrix `.npz` files to calculate electron dynamics.
