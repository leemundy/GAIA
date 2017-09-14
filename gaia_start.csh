#
#  this is a startup file for the GAIA project
#  for (t)csh users
#  Add the following to you ~/.cshrc file:
#
#       source ~/GAIA/gaia_start.csh
#
#  (or where ever the GAIA directory lives)



#  astroload comments:
#  git:
#     on the SL63 you will need a newer git to push. cloning seems to work with the old one.
#  python3:
#     we need jupyter, and perhaps astropy and a few others
#  ds9:
#     of course ds9


if (-d /astromake) then
    source /astromake/astromake_start.csh
    astroload ds9
    astroload git
    astroload -v anaconda3 python
    # R and matlab are from the default environment
else
    echo Warning, there is no /astromake here, that is odd. Assuming your environment is ok.
endif


# some useful (?) aliases

alias i   ipython
alias j   jupyter notebook


# we will maintain some sneaky aliases here, but this would imply 
# you can only use GAIA on our astro dept machines
# source code for foobar should be made available in the GAIA repo also!!

alias foobar  /home/lgm/bin/foobar
