# UnitCellJ

A Jupyter notebook for unit cell fitting


Unit-cell parameter fitting with statistical indices for outlier points.

This notebook shows how to conduct unit-cell fitting using the strategy presented in Holland and Redfurn (1997) using the `statsmodels` and `lmfit` packages.  

You need Miller index and two theta angle for the input.

Often unit cell fitting is an iterative process.  Peak positions obtained from peak fitting can be biased by a range of issues, including peak overlap and weak intensity.  Therefore, it is important to know which peaks are outtliers and which peaks are not.

Holland and Redfern proposed a robust statistical approach for this problem (read below).

T. J. B. Holland and S. A. T. Redfern (1997) "Unit cell refinement from powder diffraction data: the use of regression diagnostics". Mineralogical Magazine 61: 65-77.

Based on this approach, Holland and Redfern made a software package, UnitCell (see URL below).  However, the software is outdated for the recent MacOS.  Therefore, this notebook can provide a good alternative.

Major functions include:

- Unit-cell fitting  
- Outlier statistics 

(Note) 2022/03/04 Monoclinic and triclinic cells are not supported

## How to cite

S.-H. Shim (2022) Unit cell fitting in Jupyter. 
