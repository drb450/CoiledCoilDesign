#Please cite "Computational prediction of coiled-coil protein gelation dynamics and structure" by Dustin Britton

#Some files will be generated for user analysis during and after the simulation:
# "scoring" file is the main output file that will generate a list of mutations made, the zeroth mutation is the first and so on out of 1000 (or the total may be changed by the user). Following columns will be the Rosetta Score, N-temrinal electrostatic potential surface energy, C-terminal electrostatic potential surface energy, residue position, residue mutation, followed by the sequence with residues separated by columns. Finally the status to accept or reject is shown next to the integer probability as to whether the mutation would be accepted in the final column.
# "testfile" generates the mutation mutation made out of 1000, Rosetta score, N-temrinal electrostatic potential surface energy, C-terminal electrostatic potential surface energy of the latest mutation
# "testfile4" geneartes the status of mutation as "nope" or "huh" as well as 0 or 1 if accept or reject respecitvley. Finally, the randIndex is displayed to show the integer percent probability of accepting the mutation
# "testfile3" generates the Rosetta Score and which pose was selected between temp_input_0001.pdb-temp_input_0005.pdb of the latest mutation
# "testfile5" generates a list of the 5 Rosetta Scores generated as a result of the relax script through Rosetta of the latest mutation

# 1) The script makes a random mutation in the residuelist numbers that you provide and based on which residue list used
# 2) The script then relaxes the protein 5 times and selects the best scoring decoy for further analysis
# 3) The script then runs pdb2pqr and apbs software by creating scripts for electropotential calculation
# 4) The positive N-terminal energy and negative C-terminal energy (positive and negative based on typical patches of coiled-coil variants in the Montclare lab) is then calculated based on only the b,c,f residues that at the C and N terminus
# 5) The script then chooses to accept the mutation if both the electrostatic energies and Rosetta score are better or applies a unique probabilty of accept/reject criteria based on modified stat mech
# 6) The result is a multistate monte carlo simulation to optimize the coiled coil proteins based on the EEbcf and Rosetta score metrics

import os
import glob
import sys
import random
from shutil import copyfile
import subprocess
import math
import time


prev_score = -658  # energy score of your input variant
# GoalValues
Ngoal = 26678.18 #Target N terminal electrostatic potential energy score on the surface
Cgoal = 79814.70 #Target C terminal electrostatic potential energy score on the surface
prev_Ndif = abs(28487 - Ngoal)  # Replace values here with input N terminal bcf electrostatic potential energy score
prev_Cdif = abs(73823 - Cgoal)  # Replace values here with input C terminal bcf electrostatic potential energy score
scoreList = []
residueListNumber1 = [20,23,27,30,34,37,41,44,51] #a and d residues
residueListNumber2 = [24,26,31,33,38,40,45,47,52,54,] #e and g residues
indexRange = list(range(100))
residueList1 = ['V', 'I', 'M', 'T', 'Q', 'L']  # as listed to be favorable by woolfson
residueList2 = ['A', 'E', 'K', 'Q', 'N', 'T', 'D',
                'S']  # favored for b,c,e,f,g residues and added N,T because similar to Q and D S

zeta=3.93*10**-5 #RS factor
alpha=1.96*10**-4 #Cdif factor
beta=1.31*10**-4 #Ndif factor
for mut in range(1000):

    residueFile = open("mut.resfile", 'wt')
    ListChoice = random.choice([0, 1])
    if ListChoice == 0:
        residueListNumber = residueListNumber1
        residueList = residueList1
    else:
        residueListNumber = residueListNumber2
        residueList = residueList2
    resMut = random.choice(residueListNumber)
    residue = random.choice(residueList)
    residueFile.write("NATAA \nstart \n\n" + str(resMut) + ' A PIKAA ' + str(residue) + '\n')

    residueFile.close()

    prev_resMut = resMut
    proc1 = subprocess.call([
                                '/share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease',
                                '-parser:protocol symandrelax.xml', '-s input.pdb', '-nstruct 5', '-out:prefix temp_'])
    time.sleep(2)
    bestInput = ''
    copyfile('temp_score.sc', str(mut) + '_score.sc')
    scoreFile = open('temp_score.sc', 'r')
    lines = scoreFile.readlines()
    scoreFile.close()
    scores = []
    for i in range(len(lines)):
        line = lines[i].split()
        if i >= 2:
            scores.append(float(line[1]))
    scoreVal = min(scores)
    testy = open('testfile5', 'wt')
    testy.write(str(scores))
    testy.close()
    for i in range(len(lines)):
        line = lines[i].split()
        if i >= 2:
            if float(line[1]) == scoreVal:
                bestInput = line[21] + '.pdb'
    testy = open('testfile3', 'wt')
    testy.write(str(scoreVal) + ' ' + bestInput)
    testy.close()
    # Create an input file with a single chain in case of success
    fullprotein = open(bestInput, 'r')
    chainprotein = open('output_chain.pdb', 'wt')
    temporaryFile = open('tempest.pdb', 'wt')
    tercount = 0
    for line in fullprotein:
        if line.count('TER') == 0:
            chainprotein.write(line)
        else:
            break
    fullprotein.close()
    fullprotein = open(bestInput, 'r')
    for line in fullprotein:
        if line.count('ATOM') == 1:
            temporaryFile.write(line)
    fullprotein.close()
    chainprotein.close()
    temporaryFile.close()
    os.remove(bestInput)
    os.rename('tempest.pdb', bestInput)
    bestChain = 'output_chain.pdb'

    proc2 = subprocess.call(['/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30', '--apbs-input=output_pqr.in', '--with-ph=8',
                             '--titration-state-method=propka', '--ff=AMBER', bestInput, 'output_pqr'])
    for file in glob.glob("*.err"):
        check_errFile=file
    check_err = open(check_errFile,'r')
    if "TypeError: '<' not supported between instances of" in check_err.read().splitlines():
        proc2 = subprocess.call(['/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30', '--apbs-input=output_pqr.in', '--with-ph=8',
                                '--titration-state-method=propka', '--ff=AMBER', '--noopts', bestInput, 'output_pqr'])
    check_err.close()
    # Modify PQR output for full residue calculation
    time.sleep(2)
    dotIn = open('output_pqr.in', 'r')
    dotIn2 = open('output.in', 'w')
    for line in dotIn:
        line = line.replace('total', 'comps')
        dotIn2.write(line)
    dotIn.close()
    dotIn2.close()
    proc3 = subprocess.call(['apbs', '--output-file=output_apbs', 'output.in'])

    fpqr = open('output_pqr', 'r')
    # Initialize lists
    lines = fpqr.readlines()

    atomNum = []
    resName = []
    resNum = []

    for i in range(len(lines)):
        line = lines[i].split()
        if line[0] == 'ATOM':
            resName.append(line[3])
            atomNum.append(line[1])
            resNum.append(line[4])
    fpqr.close()

    # Read abps file
    fapbs = open('output_apbs', 'r')
    lines = fapbs.readlines()
    fapbs.close()
    polarEE = []
    elecE = []
    for i in range(len(lines)):
        line = lines[i].split()
        if i >= 2:
            if line[0] == 'atom':
                elecE.append(float(line[2]))
    polarEE = elecE[len(elecE) / 2:len(elecE)]
    apolarEE = elecE[0:len(elecE) / 2]
    NegEng = sum(polarEE)
    FIRST = '17'
    LAST = '55'
    helix_count = 5
    total_count = 0

    # b,c,f electropotential energy calculation
    # Positivie Energy
    seq = 'V '
    resAdd = ''
    posEng = 0
    negEng = 0
    stop_spot = 11
    name_temp = FIRST
    nbcfSeq = ''
    nbcfFile = open('nbcfFile', 'w')
    for i in range(0, len(polarEE) / 5):
        if resNum[i] != name_temp:
            helix_count = helix_count + 1
            total_count = total_count + 1
            name_temp = resNum[i]
            if resName[i] == 'GLY':
                resAdd = 'G '
            if resName[i] == 'ALA':
                resAdd = 'A '
            if resName[i] == 'SER':
                resAdd = 'S '
            if resName[i] == 'THR':
                resAdd = 'T '
            if resName[i] == 'CYS':
                resAdd = 'C '
            if resName[i] == 'VAL':
                resAdd = 'V '
            if resName[i] == 'LEU':
                resAdd = 'L '
            if resName[i] == 'ILE':
                resAdd = 'I '
            if resName[i] == 'MET':
                resAdd = 'M '
            if resName[i] == 'PRO':
                resAdd = 'P '
            if resName[i] == 'PHE':
                resAdd = 'F '
            if resName[i] == 'TYR':
                resAdd = 'Y '
            if resName[i] == 'TRP':
                resAdd = 'W '
            if resName[i] == 'ASP':
                resAdd = 'D '
            if resName[i] == 'GLU':
                resAdd = 'E '
            if resName[i] == 'ASN':
                resAdd = 'N '
            if resName[i] == 'GLN':
                resAdd = 'Q '
            if resName[i] == 'HIS':
                resAdd = 'H '
            if resName[i] == 'LYS':
                resAdd = 'K '
            if resName[i] == 'ARG':
                resAdd = 'R '
            seq = seq + resAdd

            if helix_count == 8:
                helix_count = 1
        if total_count <= stop_spot:
            if helix_count == 2 or helix_count == 3 or helix_count == 6:
                posEng = posEng + polarEE[i]
                nbcfSeq = nbcfSeq + ' ' + name_temp
    nbcfFile.write(nbcfSeq + '\n')

    # Negative Energy
    helix_count = 8
    total_count = 0
    name_temp = LAST
    revpolarEE = polarEE[::-1]
    revResNum = resNum[::-1]
    name_temp = LAST
    stop_spot=26
    for i in range(0, len(revpolarEE) / 5):
        if revResNum[i] != name_temp:
            helix_count = helix_count - 1
            name_temp = revResNum[i]
            total_count = total_count + 1
            if helix_count == 0:
                helix_count = 7
        if total_count <= stop_spot:
            if helix_count == 2 or helix_count == 3 or helix_count == 6:
                negEng = negEng + revpolarEE[i]
                nbcfSeq = nbcfSeq + ' ' + name_temp
    nbcfFile.write(nbcfSeq + '\n')
    nbcfFile.close()
    PosEng = (posEng) * 5
    NegEng = negEng * 5

    numTemp = FIRST

    # END OF READ PQR READ PY FILE
    testFile = open("testfile", 'wt')
    testFile.write(str(mut) + "/" + str('1000') + " " + str(scoreVal) + " " + str(PosEng) + " " + str(NegEng))
    testFile.close()
    tracking = open('scoring', 'a+')
    tracking.write(
        str(mut) + "/" + str('1000') + " " + str(scoreVal) + " " + str(PosEng) + ' ' + str(NegEng) + ' ' + str(
            resMut) + ' ' + str(residue) + ' ' + seq)
    tracking.close()
    huh = 'nope'
    scoreList.append(str(mut) + " " + str(scoreVal) + " " + str(PosEng) + " " + str(NegEng) + " " + str(resMut) + " ")

    # DifValues

    Ndif = abs(PosEng - Ngoal)
    Cdif = abs(NegEng - Cgoal)

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
            probablity = 0.0
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
            probablity = 0.0
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
            probablity = 0.0
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
            probablity1 = 0.0

        delE2 = scoreVal - prev_score
        expFact = delE2 / (prev_score * zeta)
        try:
            probability2 = math.exp(expFact)
        except OverflowError:
            probablity2 = 0.0

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
            probablity1 = 0.0

        delE2 = scoreVal - prev_score
        expFact = delE2 / (prev_score * zeta)
        try:
            probability2 = math.exp(expFact)
        except OverflowError:
            probablity2 = 0.0

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
            probablity1 = 0.0

        delE2 = Cdif - prev_Cdif
        expFact = -delE2 / (prev_Cdif * beta)
        try:
            probability2 = math.exp(expFact)
        except OverflowError:
            probablity2 = 0.0

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
            probablity1 = 0.0

        delE2 = scoreVal - prev_score
        expFact = delE2 / (prev_score * zeta)
        try:
            probability2 = math.exp(expFact)
        except OverflowError:
            probablity2 = 0.0

        delE3 = Cdif - prev_Cdif
        expFact = -delE3 / (prev_Cdif * alpha)
        try:
            probability3 = math.exp(expFact)
        except OverflowError:
            probablity3 = 0.0

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
        tracking.write(' '+'accept '+str(probPick)+ '\n')
    else:
        tracking.write(' '+'reject '+str(probPick)+ '\n')
    tracking.close()
    testy = open('testfile4', 'wt')
    testy.write(huh + ' probability: ' + str(probPick) + ' was it a 1 or 0: ' + str(accept) + ' the randIndex was: ' + str(
        randIndex))
    testy.close()
