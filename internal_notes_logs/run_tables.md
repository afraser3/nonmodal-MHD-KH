# Tables of various run parameters
Listing here the various parameters relevant to the .cfg files found in python/config_files/

Following Jeff's suggestion of letting lettered runs (run_A.cfg, run_B.cfg, etc) denote startup/exploratory runs, 
while numbered runs will denote science/publication runs. 

In the second column, "start", "stop", and "num" are describing the scan in MA that I did: 
it starts at "start", ends at "stop", and there were "num" points in the scan. 

## Lettered runs
| Letter | MA (start, stop, num) | Re     | Pm   | kx   |  k   | Nz   | mu  | Lz/pi | Comments |
| ------ | --------------------- | ---    | ---  | ---  | ---  | ---  | --- | ----- | -------- |
| A      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 200  | 256  | 0     | 10    | Converged|
| B      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 100  | 256  | 0     | 10    |          |
| C      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 200  | 512  | 0     | 10    |          |
| D      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 400  | 512  | 0     | 10    |          |
| E      | (0.05, 2.0, 40)       | 50     | 1    | 0.4  | 200  | 256  | 0     | 10    |          |
| F      | (0.05, 2.0, 40)       | 50     | 1    | 0.6  | 200  | 256  | 0     | 10    |          |
| G      | (0.05, 2.0, 40)       | 50     | 1    | 0.8  | 200  | 256  | 0     | 10    |          |
| H      | (0.05, 2.0, 40)       | 50     | 0.5  | 0.2  | 200  | 256  | 0     | 10    |          |
| I      | (0.05, 2.0, 40)       | 50     | 0.1  | 0.2  | 200  | 256  | 0     | 10    |          |
| J      | (0.05, 2.0, 40)       | 50     | 0.05 | 0.2  | 200  | 256  | 0     | 10    |          |
| K      | (0.05, 2.0, 40)       | 50     | 0.01 | 0.2  | 200  | 256  | 0     | 10    |          |
| L      | (0.05, 2.0, 40)       | 50     | 10   | 0.2  | 200  | 256  | 0     | 10    | Convergence sketchy |
| M      | (0.05, 2.0, 40)       | 50     | 100  | 0.2  | 200  | 256  | 0     | 10    | Not converged |
| N      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 400  | 1024 | 0     | 10    |          |
| O      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 1024 | 0     | 10    | Hugely different from N at highest MA |
| P      | (0.05, 2.0, 40)       | 5      | 10   | 0.2  | 200  | 256  | 0     | 10    |          |
| Q      | (0.05, 2.0, 40)       | 0.5    | 100  | 0.2  | 200  | 256  | 0     | 10    |          |
| R      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 2048 | 0     | 10    | Similar to O, but still weird -- very jumpy wrt MA and not at same MA as O|
| S      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1000 | 1024 | 0     | 10    | Still different and weird |
| T      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 2048 | 1024 | 0     | 10    | Strange errors: x-th leading minor of the array is not positive-definite in scipy.linalg.cholesky call |
| U      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1024 | 512  | 0     | 10    | Same error |
| V      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1000 | 512  | 0     | 10    | Same error |
| W      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 600  | 512  | 0     | 10    | No error; "jumpiness" persists |
| X      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 512  | 0     | 10    | No error; "jumpiness" persists |
| Y      | (0.05, 4.0, 80)       | 50     | 1    | 0.2  | 200  | 256  | 0     | 10    | Nothing suspicious |
| Z      | (0.05, 4.0, 80)       | 50     | 1    | 0.2  | 400  | 256  | 0     | 10    | Nothing suspicious |
| AA     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1600 | 2048 | 0     | 10    |          |
| AB     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1600 | 1024 | 0     | 10    |          |
| AC     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 2000 | 2048 | 0     | 10    |          |
| AD     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 1024 | 0.25  | 10    | Oops, should've done complex mu |
| AE     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 1024 | 0.25j | 10    | This is an improvement over O |
| AF     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 1024 | 0.15j | 10    | Maybe slightly better than AF |
| AG     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 400  | 512  | 0.15j | 10    | Best-looking one yet |
| AH     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 512  | 0.15j | 10    | Not quite as good |
| AI     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 512  | 1.0j  | 10    | Pretty bad |
| AJ     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 512  | 10j   | 10    | Very very bad |
| AK     | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 512  | -1.0j | 10    | A new level of bad |
| AL     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 400  | 512  | 0.15j | 10    | Not converged |
| AM     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 400  | 1024 | 0.15j | 10    | Probably not converged |
| AN     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 800  | 1024 | 0.15j | 10    | Probably not converged |
| AO     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 400  | 2048 | 0.1j  | 10    |  |
| AP     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 800  | 2048 | 0.1j  | 10    |  |
| AQ     | (0.05, 4.0, 80)       | 2500   | 1    | 0.2  | 1600 | 2048 | 0.1j  | 10    |  |