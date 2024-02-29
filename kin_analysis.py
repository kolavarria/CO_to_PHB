# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:39:04 2022

@author: josch
"""
from equilibrator_api import Q_ #for creating quantity variables
import equilibrator_custom_functions_my as eq #custom functions created during analysis
import numpy as np #numerical package
import math
import matplotlib.pyplot as plt #plotting package
import matplotlib.ticker as mtick

import custom_plot_functions as cpf

# if more space is required
from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))
display(HTML("<style>div.output_scroll {height: 50em; }</style>"))

# set defaults for plotting - grid on, linewidth = 2
plt.rc( 'axes', grid=True  )
plt.rc( 'figure', figsize = (7,4), dpi=96)
plt.rc( 'axes', linewidth=1 )
plt.rc( 'lines', linewidth=2 )

#E. coli physiological conditions as a dictionary
E_coli_con = {
            "p_h": str(Q_("7.5")),
            "p_mg": str(Q_(10)),
            "ionic_strength": str(Q_('200 mM')),
            "temperature": str(Q_(37+273.15,"K")),
            "kcat_source": "fwd",
}

#Generating Compound Settings with the default eQuilibrator bounds
cs_default_bounds = eq.obtain_compound_settings("metabolite_reference_table", custom_bounds = False)
#Lowering the CoA lower concentration bound to 1 uM, 'freeing' CoA concentration for optimization at lower levels
cs_free_CoA = eq.change_bounds(cs_default_bounds, [('CoA',Q_(1e-6,'M'), Q_(1e-3,'M'))])
#Reducing AcCoA upper bound in order to represent conditions of lower size of CoA moiety of approximately 1 mM
cs_reduced_AcCoA_CoA = eq.change_bounds(cs_free_CoA, [('AcCoA',Q_(1e-6,'M'), Q_(1e-3,'M'))])
#NAD and NADH are released to bounds of 0.1 mM and 1 mM
cs_free_NAD = eq.change_bounds(cs_free_CoA, [('NAD',Q_(1e-4,'M'),Q_(1e-3,'M')),('NADH',Q_(1e-4,'M'),Q_(1e-3,'M'))])

CDW__Cmol_ratio = Q_(24.56,'g/mol')

Volume_per_CDW = Q_(1.9, 'ul/mg').to('L/g')
Y_ATP = Q_(2,'mol/mol') #mol ATP/mol pathway net reaction
mATP_perCDW = Q_(3.2e-3,'mol/gram/hour')

flux_ATP = mATP_perCDW/Y_ATP/Volume_per_CDW
flux_ATP = flux_ATP.to('M/s')
print(flux_ATP)

##
(Keq,_) = eq.check_parameters("NOGEMP",cs_default_bounds,E_coli_con,flux_ATP, tolerance = 0.05)