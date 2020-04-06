import sys

na = 10
nb = 9
spin = 2

print "<input>"
print ""

print"  <cirun>"
print ""

print "    <general>"
print "      <spin>", spin, "</spin>"
print "      <electron na='", na, "' nb='", nb, "' />"
print "    </general>"

print ""

print "    <block>"
print "      <space ne='", 4, "' imo='", 1, "' fmo='", 2, "' />"
print "      <space ne='", 1, "' imo='", 3, "' fmo='", 3, "' />"
print "      <space ne='", 12, "' imo='", 4, "' fmo='", 9, "' />"
print "      <space ne='", 2, "' imo='", 10, "' fmo='", 11, "' />"
print "    </block>"

print ""
print "  </cirun>"
print ""

for i in range(12,200):

  print"  <cirun>"

  print "    <general>"
  print "      <spin>", spin, "</spin>"
  print "      <electron na='", na, "' nb='", nb, "' />"
  print "    </general>"

  print ""

  print "    <block>"
  print "      <space ne='", 4, "' imo='", 1, "' fmo='", 2, "' />"
  print "      <space ne='", 2, "' imo='", 3, "' fmo='", 3, "' />"
  print "      <space ne='", 12, "' imo='", 4, "' fmo='", 9, "' />"
  print "      <space ne='", 1, "' imo='", i, "' fmo='", i, "' />"
  print "    </block>"

  print "  </cirun>"
  print ""

print "</input>"

