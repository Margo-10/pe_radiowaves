# pe_radiowaves

The main.c contains code based on the finite difference method: the Crank-Nicolson scheme and the Thomas algorithm are applied. The input data is a file parameters.txt with wave, grid and environment parameters given in the SI system, for each parameter there is a comment inside the file. The output is a binary file output.bin containing the values ​​of the field modulus on the entire computational grid.

In the interpretation.ipynb, data from output.bin file is analyzed, calculations are carried out to determine how the relative permittivity depends on different parameters of the dust environment, and graphs are created based on the findings.
