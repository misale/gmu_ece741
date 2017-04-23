def sir_second(dr,n) :
   sir = 1/( (dr-1)**(-n) + (dr+1)**(-n) + ((2*dr)-1)**(-n) + ((2*dr)+1)**(-n) )
   sirdb = 10 * math.log10(sir)
   return sirdb
