# Tables of various run parameters
Listing here the various parameters relevant to the .cfg files found in python/config_files/

Following Jeff's suggestion of letting lettered runs (run_A.cfg, run_B.cfg, etc) denote startup/exploratory runs, 
while numbered runs will denote science/publication runs. 

In the second column, "start", "stop", and "num" are describing the scan in MA that I did: 
it starts at "start", ends at "stop", and there were "num" points in the scan. 

## Lettered runs
| Letter | MA (start, stop, num) | Re     | Pm   | kx   |  k   | Nz   | Lz/pi | Comments |
| ------ | --------------------- | ---    | ---  | ---  | ---  | ---  | ----- | -------- |
| A      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 200  | 256  | 10    | Converged|
| B      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 100  | 256  | 10    |          |
| C      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 200  | 512  | 10    |          |
| D      | (0.05, 2.0, 40)       | 50     | 1    | 0.2  | 400  | 512  | 10    |          |
| E      | (0.05, 2.0, 40)       | 50     | 1    | 0.4  | 200  | 256  | 10    |          |
| F      | (0.05, 2.0, 40)       | 50     | 1    | 0.6  | 200  | 256  | 10    |          |
| G      | (0.05, 2.0, 40)       | 50     | 1    | 0.8  | 200  | 256  | 10    |          |
| H      | (0.05, 2.0, 40)       | 50     | 0.5  | 0.2  | 200  | 256  | 10    |          |
| I      | (0.05, 2.0, 40)       | 50     | 0.1  | 0.2  | 200  | 256  | 10    |          |
| J      | (0.05, 2.0, 40)       | 50     | 0.05 | 0.2  | 200  | 256  | 10    |          |
| K      | (0.05, 2.0, 40)       | 50     | 0.01 | 0.2  | 200  | 256  | 10    |          |
| L      | (0.05, 2.0, 40)       | 50     | 10   | 0.2  | 200  | 256  | 10    | Convergence sketchy |
| M      | (0.05, 2.0, 40)       | 50     | 100  | 0.2  | 200  | 256  | 10    | Not converged |
| N      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 400  | 1024 | 10    |          |
| O      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 1024 | 10    | Hugely different from N at highest MA |
| P      | (0.05, 2.0, 40)       | 5      | 10   | 0.2  | 200  | 256  | 10    |          |
| Q      | (0.05, 2.0, 40)       | 0.5    | 100  | 0.2  | 200  | 256  | 10    |          |
| R      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 800  | 2048 | 10    |          |
| S      | (0.05, 4.0, 80)       | 500    | 1    | 0.2  | 1000 | 1024 | 10    |          |