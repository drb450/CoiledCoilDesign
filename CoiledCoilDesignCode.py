# This script implements a Monte Carlo simulation for optimizing coiled-coil proteins.
# Please cite "Computational prediction of coiled-coil protein gelation dynamics and structure" by Dustin Britton.

# Necessary imports
import os
import glob
import sys
import random
from shutil import copyfile
import subprocess
import math
import time

# Initial energy score and goals
prev_score = -658  # Initial Rosetta score of the input variant
Ngoal = 26678.18  # Target N-terminal electrostatic potential energy
Cgoal = 79814.70  # Target C-terminal electrostatic potential energy

# Initial electrostatic energy differences
prev_Ndif = abs(28487 - Ngoal)
prev_Cdif = abs(73823 - Cgoal)

# Residue lists and probabilities for mutation
scoreList = []
residueListNumber1 = [20, 23, 27, 30, 34, 37, 41, 44, 51]  # a and d residues
residueListNumber2 = [24, 26, 31, 33, 38, 40, 45, 47, 52, 54]  # e and g residues
residueList1 = ['V', 'I', 'M', 'T', 'Q', 'L']  # Favorable residues by Woolfson
residueList2 = ['A', 'E', 'K', 'Q', 'N', 'T', 'D', 'S']  # Favored for b, c, e, f, g residues

# Scaling factors for accept/reject probability
zeta = 3.93 * 10**-5
alpha = 1.96 * 10**-4
beta = 1.31 * 10**-4

# Monte Carlo simulation loop
for mut in range(1000):
    # Create mutation resfile
    residueFile = open("mut.resfile", 'wt')
    ListChoice = random.choice([0, 1])
    residueListNumber = residueListNumber1 if ListChoice == 0 else residueListNumber2
    residueList = residueList1 if ListChoice == 0 else residueList2
    resMut = random.choice(residueListNumber)
    residue = random.choice(residueList)
    residueFile.write(f"NATAA \nstart \n\n{resMut} A PIKAA {residue}\n")
    residueFile.close()

    # Run Rosetta for protein relaxation and scoring
    subprocess.call([
        '/share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease',
        '-parser:protocol', 'symandrelax.xml', '-s', 'input.pdb', '-nstruct', '5', '-out:prefix', 'temp_'
    ])
    time.sleep(2)

    # Analyze scores and determine the best scoring structure
    scores = []
    bestInput = ''
    with open('temp_score.sc', 'r') as scoreFile:
        lines = scoreFile.readlines()
    for line in lines[2:]:  # Skip header lines
        scores.append(float(line.split()[1]))
    scoreVal = min(scores)

    with open('testfile5', 'wt') as testy:
        testy.write(str(scores))

    for line in lines[2:]:
        if float(line.split()[1]) == scoreVal:
            bestInput = line.split()[21] + '.pdb'

    with open('testfile3', 'wt') as testy:
        testy.write(f"{scoreVal} {bestInput}")

    # Process the best scoring structure
    with open(bestInput, 'r') as fullprotein, open('output_chain.pdb', 'wt') as chainprotein:
        for line in fullprotein:
            if 'TER' not in line:
                chainprotein.write(line)
            else:
                break

    subprocess.call([
        '/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30',
        '--apbs-input=output_pqr.in', '--with-ph=8', '--titration-state-method=propka', '--ff=AMBER', bestInput,
        'output_pqr'
    ])

    # Handle potential errors from PDB2PQR
    for err_file in glob.glob("*.err"):
        with open(err_file, 'r') as check_err:
            if "TypeError: '<' not supported between instances of" in check_err.read():
                subprocess.call([
                    '/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30',
                    '--apbs-input=output_pqr.in', '--with-ph=8', '--titration-state-method=propka', '--ff=AMBER',
                    '--noopts', bestInput, 'output_pqr'
                ])

    # Modify input for APBS and calculate electrostatic energies
    with open('output_pqr.in', 'r') as dotIn, open('output.in', 'w') as dotIn2:
        for line in dotIn:
            dotIn2.write(line.replace('total', 'comps'))

    subprocess.call(['apbs', '--output-file=output_apbs', 'output.in'])

    # Read energies and calculate properties
    with open('output_apbs', 'r') as fapbs:
        lines = fapbs.readlines()
    elecE = [float(line.split()[2]) for line in lines if line.startswith('atom')]
    polarEE = elecE[len(elecE)//2:]
    apolarEE = elecE[:len(elecE)//2]

    posEng = sum(polarEE[:5])
    negEng = sum(polarEE[-5:])
    
    # Record results
    with open('testfile', 'wt') as testFile:
        testFile.write(f"{mut}/1000 {scoreVal} {posEng} {negEng}")

    # Scoring and acceptance condition
    Ndif = abs(posEng - Ngoal)
    Cdif = abs(negEng - Cgoal)
    # accept Condition
	if scoreVal <= prev_score and Ndif <= prev_Ndif and Cdif <= prev_Cdif:  # changed to include PosEng score too and NegEng
		copyfile(bestInput, str(mut) + '_mut1_output.pdb')
		probPick = 1
		accept = 1
		randIndex = 1
		prev_score = scoreVal
		prev_Ndif = Ndif
		prev_Cdif = Cdif
		os.rename(bestChain, 'input.pdb')
		os.remove(bestInput)
		huh = 'mutated'
	# R condition - score is bad
	elif scoreVal > prev_score and Ndif <= prev_Ndif and Cdif <= prev_Cdif:
		delE = scoreVal - prev_score
		expFact = delE / (prev_score * zeta)
		try:
			probability = math.exp(expFact)
		except OverflowError:
			probability = 0.0
		probPick = int(round(probability * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			huh = 'mutated'
			copyfile(bestInput, str(mut) + '_mutR_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestChain)
	# P Condition - N term is bad
	elif scoreVal <= prev_score and Ndif > prev_Ndif and Cdif <= prev_Cdif:
		delE = Ndif - prev_Ndif
		expFact = -delE / (prev_Ndif * beta)
		try:
			probability = math.exp(expFact)
		except OverflowError:
			probability = 0.0
		probPick = int(round(probability * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			huh = 'mutated'
			copyfile(bestInput, str(mut) + '_mutC_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)
	# N condition - C term is bad
	elif scoreVal <= prev_score and Ndif <= prev_Ndif and Cdif > prev_Cdif:
		delE = Cdif - prev_Cdif
		expFact = -delE / (prev_Cdif * alpha)
		try:
			probability = math.exp(expFact)
		except OverflowError:
			probability = 0.0
		probPick = int(round(probability * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			huh = 'mutated'
			copyfile(bestInput, str(mut) + '_mutT_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)
	# R condition x P condition - score and N term is bad
	elif scoreVal > prev_score and Ndif > prev_Ndif and Cdif <= prev_Cdif:
		delE1 = Ndif - prev_Ndif
		expFact = -delE1 / (prev_Ndif * beta)
		try:
			probability1 = math.exp(expFact)
		except OverflowError:
			probability1 = 0.0

		delE2 = scoreVal - prev_score
		expFact = delE2 / (prev_score * zeta)
		try:
			probability2 = math.exp(expFact)
		except OverflowError:
			probability2 = 0.0

		probPick = int(round(probability1 * probability2 * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			copyfile(bestInput, str(mut) + '_mutRC_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
			huh = 'mutated'
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)
	# R condition x N condition
	elif scoreVal > prev_score and Ndif <= prev_Ndif and Cdif > prev_Cdif:
		delE1 = Cdif - prev_Cdif
		expFact = -delE1 / (prev_Cdif * alpha)
		try:
			probability1 = math.exp(expFact)
		except OverflowError:
			probability1 = 0.0

		delE2 = scoreVal - prev_score
		expFact = delE2 / (prev_score * zeta)
		try:
			probability2 = math.exp(expFact)
		except OverflowError:
			probability2 = 0.0

		probPick = int(round(probability1 * probability2 * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			copyfile(bestInput, str(mut) + '_mutRT_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
			huh = 'mutated'
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)
	# P condition x N condition
	elif scoreVal <= prev_score and Ndif > prev_Ndif and Cdif > prev_Cdif:
		delE1 = Ndif - prev_Ndif
		expFact = -delE1 / (prev_Ndif * alpha)
		try:
			probability1 = math.exp(expFact)
		except OverflowError:
			probability1 = 0.0

		delE2 = Cdif - prev_Cdif
		expFact = -delE2 / (prev_Cdif * beta)
		try:
			probability2 = math.exp(expFact)
		except OverflowError:
			probability2 = 0.0

		probPick = int(round(probability1 * probability2 * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			copyfile(bestInput, str(mut) + '_mutCT_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
			huh = 'mutated'
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)
	# R condition x P condition x N condition
	elif scoreVal > prev_score and Ndif > prev_Ndif and Cdif > prev_Cdif:
		delE1 = Ndif - prev_Ndif
		expFact = -delE1 / (prev_Ndif * beta)
		try:
			probability1 = math.exp(expFact)
		except OverflowError:
			probability1 = 0.0

		delE2 = scoreVal - prev_score
		expFact = delE2 / (prev_score * zeta)
		try:
			probability2 = math.exp(expFact)
		except OverflowError:
			probability2 = 0.0

		delE3 = Cdif - prev_Cdif
		expFact = -delE3 / (prev_Cdif * alpha)
		try:
			probability3 = math.exp(expFact)
		except OverflowError:
			probability3 = 0.0

		probPick = int(round(probability1 * probability2 * probability3 * 100, 0))
		pickList = [0] * 100
		for j in range(len(pickList)):
			if j < probPick:
				pickList[j] = 1
		randIndex = random.choice(indexRange)
		accept = pickList[randIndex]
		if accept == 1:
			copyfile(bestInput, str(mut) + '_mutRCT_output.pdb')
			prev_score = scoreVal
			prev_Ndif = Ndif
			prev_Cdif = Cdif
			os.rename(bestChain, 'input.pdb')
			os.remove(bestInput)
			huh = 'mutated'
		else:
			copyfile(bestInput, str(mut) + '_output.pdb')
			os.remove(bestInput)

	for fname in os.listdir('./'):
		if fname.startswith("output") or fname.startswith("temp_"):
			os.remove(fname)
	tracking = open('scoring', 'a+')
	if huh == 'mutated':
		tracking.write(' ' + 'accept ' + str(probPick) + '\n')
	else:
		tracking.write(' ' + 'reject ' + str(probPick) + '\n')
	tracking.close()
	testy = open('testfile4', 'wt')
	testy.write(huh + ' probability: ' + str(probPick) + ' was it a 1 or 0: ' + str(accept) + ' the randIndex was: ' + str(randIndex))
	testy.close()
