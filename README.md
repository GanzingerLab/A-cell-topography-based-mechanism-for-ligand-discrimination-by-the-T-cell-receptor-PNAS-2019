# A-cell-topography-based-mechanism-for-ligand-discrimination-by-the-T-cell-receptor-PNAS-2019
This is the simulation code published in "A cell topography-based mechanism for ligand discrimination by the T cell receptor; Ricardo A. Fernandes, Kristina A. Ganzinger, Justin C. Tzou, Peter Jönsson, Steven F. Lee, Matthieu Palayret, Ana Mafalda Santos, Alexander R. Carr, Aleks Ponjavic, Veronica T. Chang, Charlotte Macleod, B. Christoffer Lagerholm, Alan E. Lindsay, Omer Dushek, Andreas Tilevik, Simon J. Davis, David Klenerman; Proceedings of the National Academy of Sciences Jul 2019, 116 (28) 14002-14010; DOI: 10.1073/pnas.1817255116"

*Please cite the paper if you use it.*

The article can be downloaded (open access) under https://www.pnas.org/content/116/28/14002

Briefly, the code (published model moving boundaries) is structured by paper figures: all figures from the paper can be regenerated using the "GetPm" function. This functions loads solutions from PDEs that are calculated by 'MFPTGrowingDisk.m' code (which is computationally expensuive, so running this on a cluster is sensible) - this fucntion uses 'MFPTLoop.m' to calculate the solutions for many runs.

The modular nature of this calculation means that some parameters, such as the rate of TCR entry, diffusion coefficients and concentrations can varied in the Pm function, but for changing other parameters, the Ps has to be re-calculated using the 'MFPTGrowingDisk.m' code (such as agonist density, k_on, k_off ). A mathematical description can be found in the the pdf "mathematical notes-publishedmodel".



