from _sw import *

def readMatrixFile(fileName, style="protein"):
	mfile = open(fileName, "r")

	matrix = {}
	if style == "protein":
		minChars = 20
	else:
		minChars = 4
	divisor = 1.0
	while 1: # skip through header
		line = mfile.readline()
		if not line:
			raise ValueError, "No matrix found in %s" % fileName
		chars = line.split()
		if chars[:2] == ["#", "Divisor:"]:
			divisor = float(chars[2])
			if divisor <= 0:
				raise ValueError("Bad divisor (%g) in matrix %s"
					% (divisor, fileName))
		if len(chars) < minChars:
			continue
		if max(map(len, chars)) > 1:
			continue

		for char in chars:
			fields = mfile.readline().split()
			if len(fields) != len(chars) + 1:
				raise ValueError, "Unexpected number of "\
					"fields on line %d of matrix in %s" % (
					chars.index(char)+1, fileName)
			if fields[0] != char:
				raise ValueError, "Line %d of %s is %s, " \
					"expected %s" % (chars.index(char)+1,
					fileName, fields[0], char)
			for i in range(len(chars)):
				matrix[(fields[0], chars[i])] = float(fields[i+1]) / divisor
		# assume '?' is an unknown residue...
		for char in chars:
			matrix[(char, '?')] = matrix[('?', char)] = 0.0
		matrix[('?', '?')] = 0.0
		# treat pyrrolysine ('O') and selenocysteine ('U') as
		# 'X' and 'C' respectively (see Trac #12828 for reasoning)
		if style == "protein":
			for char in chars + ["?"]:
				matrix[(char, 'O')] = matrix[('O', char)] = matrix[(char, 'X')]
				matrix[(char, 'U')] = matrix[('U', char)] = matrix[(char, 'C')]
			matrix[('O', 'U')] = matrix[('U', 'O')] = matrix[('X', 'C')]
			matrix[('O', 'O')] = matrix[('X', 'X')]
			matrix[('U', 'U')] = matrix[('C', 'C')]
		mfile.close()
		break
	return matrix

import chimera
import os
matrices = {}
matrixFiles = {}

for dir in chimera.pathFinder().allExistingDirs("SmithWaterman", "matrices"):
	for matrixFile in os.listdir(dir):
		if matrixFile.startswith("."):
			# some editors and file system mounters create files that start
			# with '.' that should be ignored
			continue
		if matrixFile.endswith(".matrix"):
			ftype = "protein"
			name = matrixFile[:-7]
		elif matrixFile.endswith(".dnamatrix"):
			ftype = "dna"
			name = matrixFile[:-10]
		else:
			continue
		path = os.path.join(dir, matrixFile)
		try:
			matrices[name] = readMatrixFile(path, style=ftype)
		except ValueError:
			chimera.replyobj.reportException()
		matrixFiles[name] = path
try:
	del dir, matrixFile, ftype, name
except NameError:
	chimera.replyobj.warning("No matrices found by SmithWaterman module.\n")
