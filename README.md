# FK-FGA-
Development of parquet dual fermion code for the Falicov-Kimball model.

Additional scripts supplementing the master and ladder branch.

### Calculation of chemical potential to fix occupation of mobile c electrons

- **FindMu.py**
- **CalcNewOcc.py**
- **routines.py**

If you want to investigate a system away from half-filling, you first need to find the right chemical potential that fixes the occupation of the mobile c electrons (the occupation of f electrons is fixed in the code anyhow). This is done in **FindMu.py** employing a simple bisection method where in every step the DMFT calculation **SCDMFT.out** is called.

In **CalcNewOcc.py** the resulting occupation of the DF calculation is calculated, to check how much the DMFT occupation is changed by this DF calculation.

Routines that are used in **FindMu.py** and **CalcNewOcc.py** are written in the module **routines.py**.

### Calculation of bubble term of current-current correlation

- **CalcCurrentBubble.py**
- **BubbleCurrent.py**

The calculation of the bubble term in the parquet code suffers from an incorrect frequency asymptotic due to the frequency summation in a finite frequency box. To overcome this, the bubble term can be calculated anew using **CalcCurrentBubble.py** which includes the module in **BubbleCurrent.py**. Here, the bubble term is calculated after a Fourier transform of the DF Green's function to imaginary times instead on the Matsubara frequency axis.

### Calculation of second order diagram to dual self-energy

- **CallSecondOrder.py**
- **SecondOrderFunctions.py**

Python code to calculate the second order diagram to the dual self-energy, calculated in real space and then transformed back.
