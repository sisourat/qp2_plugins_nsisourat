import sys

na = 5
nb = 4
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
print "      <space ne='", 1, "' imo='", 1, "' fmo='", 1, "' />"
print "      <space ne='", 8, "' imo='", 2, "' fmo='", 5, "' />"
print "    </block>"

print ""
print "  </cirun>"
print ""

for i in range(32,82):

  print"  <cirun>"

  print "    <general>"
  print "      <spin>", spin, "</spin>"
  print "      <electron na='", na, "' nb='", nb, "' />"
  print "    </general>"

  print ""

  print "    <block>"
  print "      <space ne='", 2, "' imo='", 1, "' fmo='", 1, "' />"
  print "      <space ne='", 6, "' imo='", 2, "' fmo='", 5, "' />"
  print "      <space ne='", 1, "' imo='", i, "' fmo='", i, "' />"
  print "    </block>"

  print "  </cirun>"
  print ""

print "</input>"

