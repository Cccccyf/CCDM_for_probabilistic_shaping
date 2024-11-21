# CCDM for Probabilistic Shaping in Python
>The whole coding is based on the C++ implementation of CCDM mentioned in https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping.git
### Code Structure
### Usage
```
from probabilisticShaping import probabilisticShaping
M = 256
shaping_factor = 0.01
symbol_num = 10000
ps = probabilisticShaping(M, shaping_factor, symbol_num)
symbols = ps.genSymbols()
```
The constellation of symbols generated is enlarged to keep same power as uniform distribution qam.
### Reference
[1]Edson Porto da Silva, Adolfo Fernandes Herbster. "OptiCommPy: Open-source Simulation of Fiber Optic Communications with Python", Journal of Open Source Software, 9(98), 6600, (2024) https://doi.org/10.21105/joss.06600
<br>
[2] https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping.git
