def sir_hex(dr,n) :
   sir = 1/( (dr-1)**(-n) + (dr+1)**(-n)  + 4*(dr**(-n)) )
   sirdb = 10 * math.log10(sir)
   print "SIR : ",sir
   print "SIRdB : ",sirdb
