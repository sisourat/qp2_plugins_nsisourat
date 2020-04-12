def compare2dets(det1alp,det1bet,det2alp,det2bet):
   ndiff = 0
   ndiff_alp = 0
   orbdiff_alp = []
   for i in range(len(det1alp)):
     if( not(det1alp[i]==det2alp[i]) ):
#       print i, det1alp[i],det2alp[i]
       orbdiff_alp.append(i)
       ndiff_alp +=1
       ndiff +=1
   ndiff_bet = 0
   orbdiff_bet = []
   for i in range(len(det1bet)):
     if( not(det1bet[i]==det2bet[i]) ):
#       print i, det1bet[i],det2bet[i]
       orbdiff_bet.append(i)
       ndiff_bet +=1
       ndiff +=1
   return ndiff, ndiff_alp, ndiff_bet, orbdiff_alp, orbdiff_bet

