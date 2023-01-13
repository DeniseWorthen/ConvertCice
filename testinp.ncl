begin

h=1.0
nlyr4=4
nlyr5=7

delh4 = h/tofloat(nlyr4)
delh5 = h/tofloat(nlyr5)

print(delh4+"   "+delh5)

c4x = fspan(delh4/2.0,h-(delh4/2.0),nlyr4)
c5x = fspan(delh5/2.0,h-(delh5/2.0),nlyr5)
print("layer axis defined for CICE4  ")
print(c4x)
print("  ")
print("layer axis defined for CICE5 ")
print(c5x)

end