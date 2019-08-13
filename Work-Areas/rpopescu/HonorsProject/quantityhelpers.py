from astropy import units as u
from astropy.units.quantity import Quantity

def isQuantity(q):
   """Test if input is an Astropy Quantity object
      Parameters: 
        q - input object to test
   """
   return type(q) == Quantity

def isFluxDensity(q):
   """Test if input is an Astropy Quantity with units of flux density
      Parameters: 
        q - input object to test
   """
   if isQuantity(q):
       return q.unit.is_equivalent("Jy")
   return False

def isLength(q):
   """Test if input is an Astropy Quantity with units of length
      Parameters: 
        q - input object to test
   """
   if isQuantity(q):
       return q.unit.is_equivalent("cm")
   return False

def isMagnitude(q):
   """Test if input is an Astropy Quantity with units of astronomical magnitude
      Parameters: 
        q - input object to test
   """
   try:
       return q.unit.is_equivalent("mag")
   except Exception:
       return False
