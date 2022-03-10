#!/usr/bin/env python
# coding: utf-8

# # Unit cell fitting

# created by S.-H. Shim, 2022/03/05

# This notebook shows how to conduct unit-cell fitting using the strategy presented in Holland and Redfurn (1997) using the `statsmodels` and `lmfit` packages.  
# 
# You need Miller index and two theta angle for the input.
# 
# Often unit cell fitting is an iterative process.  Peak positions obtained from peak fitting can be biased by a range of issues, including peak overlap and weak intensity.  Therefore, it is important to know which peak position data points are outtliers and which peak position data points are not.
# 
# Holland and Redfern proposed a statistical approach for this problem (read below).
# 
# T. J. B. Holland and S. A. T. Redfern (1997) "Unit cell refinement from powder diffraction data: the use of regression diagnostics". Mineralogical Magazine 61: 65-77.
# 
# Based on this approach, Holland and Redfern made a software package, UnitCell (see URL below).  However, the software is outdated for the recent MacOS.  
# 
# This notebook provides the same method using python modules.  Therefore, this notebook can provide a good alternative for those you want to conduct similar analysis in latest version MacOS with python.
# 
# URL: http://ccp14.cryst.bbk.ac.uk/ccp/web-mirrors/crush/astaff/holland/UnitCell.html

# ### Note
# 
# 
# At the moment, this notebook does not work for monoclic and triclinic cells.

# In[1]:


import numpy as np
import pandas as pd


# Key functions for unit cell fitting are provided in a separate python file, `ucfits.py`.

# In[2]:


from ucfits import *


# ## Wavelength
# 
# Provide wavelength of your XRD data in anstrom unit.

# In[3]:


wavelength = 0.4133


# ## Mg(OH)2 - an example for hexagonal

# You can just change the cell below.

# In[4]:


data_MgOH2 = [ # h k l twotheta
    [ 0.0, 0.0, 1.0, 5.692435],
    [ 1.0, 0.0, 0.0, 9.362543], 
    [ 0.0, 1.0, 1.0, 10.966361],
    [ 0.0, 1.0, 2.0, 14.766217],
    [ 1.0, 1.0, 0.0, 16.248651] ]


# The cell below converts your data to a pandas dataframe for input.

# In[5]:


df = data2df(data_MgOH2, wavelength)
df


# #### Unit cell fitting
# 
# For different symmetry, change the name of the function: `fit_*_cell` (* = cubic, tetragonal, hexagonal, and orthorhombic).
# 
# `cell_param`, `s_cell_param`, `cell_vol`, and `s_cell_vol` are fitted cell parameters, estimated error for the cell parameters, cell volume, and estimated cell volume.  You may use these lists if you need further calculations in this notebook.  They can be also found in the output.

# In[6]:


cell_param, s_cell_param, cell_vol, s_cell_vol, results =     fit_hexagonal_cell(df, wavelength, verbose=True)


# #### Some key statistical indices to read from the output above
# 
# In `[[Variables]]` part above, you can find fitting results for the unit-cell parameters and their estimated errors.
# 
# `[[Correlations]]` provide some useful information, too.  Of course, you want this parameter to be close to 0 - un-corrrelated.  But lower is of course better.

# In[7]:


make_output_table(results, df)


# The table above is more useful, particularly if you want to find potential outlier peak position data.  
# 
# #### twoth residue
# 
# This is simply $2\theta_{obs} - 2\theta_{cal}$.  Of course, if you have a large difference, the peak data point could be an outlier.  But the fitting is highly nonlinear and small $d$-spacing peaks tend to affect the fitting more.  In fact, in theory, small $d$-spacing peaks provide more reliablle measures of unit-cell parameters.  However, small $d$-spacing (or high $2\theta$) lines tend to be weak in intensity and tend to overlap with other peaks.  
# 
# #### hat
# 
# This index deals with influence of individual data points for the fit result. `hat = 0` means no influence.  A larger hat value means larger influence.
# 
# #### Rstudent
# 
# This index expected to be less than 2 at 95% confidence level.  If a peak shows `Rstudent > 2`, you should be suspicious about the reliability of the peak data point.
# 
# #### dfFits
# 
# This is a deletion diagnostic giving the change in the predicted values upon deletion of the data point as a multiple of the standard deviation for $Q^2$.
# 
# #### dfbetas
# 
# This provides a measure of how much the calculated value of each refined parameter would change if the rergression were rerrun without using the data point.

# URL for statsmodels on outliers influence: https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.OLSInfluence.html?highlight=olsinfluence#statsmodels.stats.outliers_influence.OLSInfluence

# ## Mg(OH)2 - an example for hexagonal

# In the example below, I commented out a peak data point to see the impact.

# In[8]:


data_MgOH2 = [
    # h k l twotheta
    [ 0.0, 0.0, 1.0, 5.692435],
    [ 1.0, 0.0, 0.0, 9.362543], 
    [ 0.0, 1.0, 1.0, 10.966361],
#    [ 0.0, 1.0, 2.0, 14.766217], # very high hat value
    [ 1.0, 1.0, 0.0, 16.248651] ]


# In[9]:


df = data2df(data_MgOH2, wavelength)
df


# In[10]:


cell_param, s_cell_param, cell_vol, s_cell_vol, results =     fit_hexagonal_cell(df, wavelength, verbose=True)


# In[11]:


make_output_table(results, df)


# Compare this result with one above with all the data ponts and check how much the unit-cell parameterr changed.
# 
# Note that some indices are not calculated, such as dfFits and dfBetas.  The reason is you have only 4 data points now for 2 parameters to fit -- insufficient number of data points.

# ## Stishovite - an example for tetragonal

# In[12]:


data_stv = [
    # h k l twotheta
    [1.0, 1.0, 0.0, 8.277043],
    [1.0, 0.0, 1.0, 10.776637],
    [1.0, 1.0, 1.0, 12.262541],
    [2.0, 1.0, 0.0, 13.083539],
    [2.0, 1.0, 1.0, 15.932403],	
    [2.0, 2.0, 0.0, 16.564225] ]


# Do not modify below.  This is to convert your data to pandas data frame for input

# In[13]:


df = data2df(data_stv, wavelength)
df


# In[14]:


cell_param, s_cell_param, cell_vol, s_cell_vol, results =     fit_tetragonal_cell(df, wavelength, verbose=True)


# In[15]:


make_output_table(results, df)


# So the peak \#4 is really concerning and should be checked.
