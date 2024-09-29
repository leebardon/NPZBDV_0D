Vm = Vmax(j)*photoTempFunc(j)

muDON  = PCm*DON/(DON + ksatDON(j))
muDOC  = PCm*DOC/(DOC + ksatDOC(j))

muNH4 = pmaxDIN*NH4/(NH4 + ksatDIN)*reminTempFunc

mu = MIN(Vm, yield(j)*muDOC, muDON + muNH4, muO)

# "growth" is the mu*B part of the B equation. 
growth = mu*B(j) 

# OM-limited growth rate, which determines actual inorganic uptake -> R_NC is the ratio of N to C for biomass (i.e. 1/5)
mu_Norg = MIN(PCm, yield(j)*muDOC, muDON, muO)

# actual inorganic uptake (a fraction of potential uptake):
uptakeNH4 = MIN(muNH4, MAX((mu - mu_Norg),0.))*X(j)*R_NC(j)


# OM-limited growth rates needed for uptake of OM calc
growth_Norg = mu_Norg*X(j)*R_NC(j)

# calculate the ratios of DON/DOC etc. uptake
Rup_NC  = MIN(R_NC(j), muDON/(muDOC + EPS))

# Here is uptake and excretion of inorganic nutrients (i.e. sinks in equations for DOC and DON)
uptakeDOC = growth/yield(j)
uptakeDON = uptakeDOC*Rup_NC

# "respDOC" and "respDON" are the remin. sources so go into DIC and NH4 pools 
respDOC = growth*(1 ./ yield(j)-1.)
respDON = uptakeDON - growth_Norg


# Above we've calculated the terms that go into the differential equations for biomass, OM, and inorganic nutrients. 
# I'm not adding these to equations here in the script because it is so different in darwin. 

