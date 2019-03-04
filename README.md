# FK-FGA-

Development of parquet dual fermion code for the Falicov-Kimball model.

**SCDMFT.cpp** contains the DMFT calculation that yields the input quantities for the DF parquet approach.

The DF parquet calculation is implemented in **DFParquet.cpp** which uses the classes contained in **DFParquetObj.cpp**, **DFConductivity.cpp** and **DFMagnetism.cpp** to do the parquet calculation and afterwards calculate the corresponding current-current and density-density correlation functions.

In **routines.cpp** there are some auxiliary routines for reading and writing binary files and calculating the dispersion relation and corrected Green's function (DMFT + DF).

The full DF parquet approach is executed in the Python script **CallParquet.py**, where first a compiled **SCDMFT.out** is called to generate the input for the iterative call of **DFParquet.out** followed by it.
