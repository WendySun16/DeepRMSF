inf = open("chargeGaff.defattr", "r")
outf = open("dictAdd.py", "w")
for line in inf:
	if not line.startswith("\t"):
		continue
	line = line.strip()
	osl, val = line.split("\t")
	if osl == "#":
		continue
	if not osl.startswith(":/amberName=") or osl.count("@") != 1:
		raise ValueError("Can't handle OSL: %s" % osl)
	resOSL, atomOSL = osl.split("@")
	amberName = eval(resOSL[12:])
	if atomOSL.startswith("/element="):
		atomKey = "chimera.Element('%s')" % atomOSL[9:]
	else:
		atomKey = '"%s"' % atomOSL.lower()
	print>>outf, '\t\t("%s", %s): %s,' % (amberName, atomKey, val)
inf.close()
outf.close()
