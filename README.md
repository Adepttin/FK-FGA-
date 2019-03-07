# FK-FGA-
Development of ladder dual fermion code for the Falicov-Kimball model.

**SCDMFT.cpp** contains the DMFT calculation that yields the input quantities for the DF parquet approach.

The DF ladder calculation is implemented in **DFphLadder.cpp** for the particle-hole ladder and in **DFppLadder.cpp** for the particle-particle ladder, which use the classes **DFphLadderObj.cpp** and **DFppLadderObj.cpp**  respectively to do the ladder calculation and afterwards calculate the corresponding current-current correlation function.

In **routines.cpp** there are some auxiliary routines for reading and writing binary files and calculating the dispersion relation and corrected Green's function (DMFT + DF).

The full DF parquet approach is executed in the Python script **LadderCaller.py**, where first a compiled **SCDMFT.out** is called to generate the input for the call of **DFphLadder.out** or **DFppLadder.out** followed by it.

