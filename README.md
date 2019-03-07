# FK-FGA-
Development of parquet dual fermion code for the Falicov-Kimball model.

Additional scripts supplementing the master and ladder branch.

### Calculation of bubble term of current-current correlation

The calculation of the bubble term in the parquet code suffers from an incorrect frequency asymptotic due to the frequency summation in a finite frequency box. To overcome this, the bubble term can be calculated anew using **CalcCurrentBubble.py** which includes the module in **BubbleCurrent.py**. Here, the bubble term is calculated after a Fourier transform of the DF Green's function to imaginary times instead on the Matsubara frequency axis.
