#
#  this is a startup file for the GAIA project
#  for (t)csh users
#  Add the following to you ~/.cshrc file:
#
#       source ~/GAIA/gaia_start.csh
#
#  (or where ever the GAIA directory lives)



if (-d /astromake) then
    source /astromake/astromake_start.csh
    astroload ds9
    astroload git
    astroload python
    # R and matlab are from the default environment
else
    echo Warning, there is no /astromake here, that is odd.
endif
