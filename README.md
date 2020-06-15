# kranka_ucm
Codes written at UC Merced by and for the Isborn research group. Most are compatible with `python3.x`
- `scripts/`: these mostly interface with the `Gaussian` program. There are some standalone codes for running TDCI calculations as well.
- `ML/`: these codes are written using `tensorflow` and `pytorch` libraries (to simulate deep, feed-forward static neural networks for time-series analysis). Inspired and modified by [`@hbhat4000`](https://github.com/hbhat4000)
- `examples/`: this directory contains example computations run using the scripts and codes stored in the above folder. `examples/CI` refers to the configuration interaction data-extraction (from `Gaussian DV` calculations) scripts (in `.ipynb` format) and the time-dependent  configuration interaction script (also in `.ipynb` format) that uses the extracted data to calculate electron dynamics.
