import sys
import numpy as np

elt = {1: 'Hydrogen', 2: 'Helium', 3: 'Lithium', 4: 'Beryllium', 5: 'Boron' , 6: 'Carbon' , 7: 'Nitrogen' , 8: 'Oxygen' , 9: 'Fluorine' , 10: 'Neon' , 11: 'Sodium' , 12: 'Magnesium' , 13: 'Aluminum' , 14: 'Silicon' , 15: 'Phosphorus' , 16: 'Sulfur' , 17: 'Chlorine' , 18: 'Argon' , 19: 'Potassium' , 20: 'Calcium' , 21: 'Scandium' , 22: 'Titanium' , 23: 'Vanadium' , 24: 'Chromium' , 25: 'Manganese' , 26: 'Iron' , 27: 'Cobalt' , 28: 'Nickel' , 29: 'Copper' , 30: 'Zinc' , 31: 'Gallium' , 32: 'Germanium' , 33: 'Arsenic' , 34: 'Selenium' , 35: 'Bromine' , 36: 'Krypton'} 

alpminH_S=0.001
alpmaxH_S=100
ngto_S=14

alpminH_P=0.001
alpmaxH_P=10.0
ngto_P=11

alpminH_D=0.001
alpmaxH_D=10.0
ngto_D=7

alpminH_F=0.001
alpmaxH_F=10.0
ngto_F=4

alpminH_G=0.001
alpmaxH_G=10.0
ngto_G=0

alpminH_H=0.001
alpmaxH_H=10.0
ngto_H=0


for e in elt:
 print elt[e]

 alpmin = alpminH_S
 alpmax = alpmaxH_S*e**2
 ngto = ngto_S
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'S 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'

 alpmin = alpminH_P
 alpmax = alpmaxH_P*e**2
 ngto = ngto_P
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'P 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'

 alpmin = alpminH_D
 alpmax = alpmaxH_D*e**2
 ngto = ngto_D
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'D 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'

 alpmin = alpminH_F
 alpmax = alpmaxH_F*e**2
 ngto = ngto_F
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'F 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'

 alpmin = alpminH_G
 alpmax = alpmaxH_G*e**2
 ngto = ngto_G
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'G 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'

 alpmin = alpminH_H
 alpmax = alpmaxH_H*e**2
 ngto = ngto_H
 if ngto == 1:
     beta = 1.0
 else:    
     beta = (alpmax/alpmin)**(1.0/(ngto-1))
 for i in range(1,ngto+1):
    print 'H 1'
    print '1 ',alpmin*beta**(i-1), ' 1.0'
 print

